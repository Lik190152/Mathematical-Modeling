from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .centerline import decompose_branches, extract_centerline_graph
from .geometry import (
    exact_edt,
    keep_largest_component,
    mask_volume_mm3,
    read_nifti_mask,
    resample_binary_mask,
)
from .lesion import detect_stenosis_candidates
from .view_selection import greedy_cover_views, select_top_local_views


TARGET_SPACING_MM = 0.4
LOCAL_TOP_K = 3
ANGLE_STEP_DEG = 5
AZIMUTH_GRID = list(range(0, 360, ANGLE_STEP_DEG))
ELEVATION_GRID = list(range(-60, 65, ANGLE_STEP_DEG))


def preprocess_mask(mask: np.ndarray, spacing: Tuple[float, float, float], target_spacing: float = TARGET_SPACING_MM):
    raw_mask = keep_largest_component(mask)
    raw_volume = mask_volume_mm3(raw_mask, spacing)
    resampled_mask, resampled_spacing, scale = resample_binary_mask(raw_mask, spacing, target_spacing)
    resampled_mask = keep_largest_component(resampled_mask)
    resampled_volume = mask_volume_mm3(resampled_mask, resampled_spacing)
    stats = {
        "raw_shape": tuple(int(v) for v in raw_mask.shape),
        "raw_spacing_mm": tuple(float(v) for v in spacing),
        "raw_volume_mm3": float(raw_volume),
        "processed_shape": tuple(int(v) for v in resampled_mask.shape),
        "processed_spacing_mm": tuple(float(v) for v in resampled_spacing),
        "processed_volume_mm3": float(resampled_volume),
        "volume_delta_pct": float(100.0 * (resampled_volume - raw_volume) / raw_volume) if raw_volume else 0.0,
    }
    return resampled_mask, resampled_spacing, stats


def sample_branch_profile(xyz_mm: np.ndarray, radii_mm: np.ndarray, step_mm: float = TARGET_SPACING_MM):
    if len(xyz_mm) == 0:
        return np.empty((0, 3), dtype=float), np.empty((0,), dtype=float), np.empty((0,), dtype=float)
    if len(xyz_mm) == 1:
        return xyz_mm.copy(), radii_mm.copy(), np.array([0.0], dtype=float)

    seg = np.linalg.norm(np.diff(xyz_mm, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(seg)])
    total_len = float(arc[-1])
    sample_arc = np.arange(0.0, total_len + step_mm * 0.5, step_mm, dtype=float)
    if sample_arc[-1] < total_len:
        sample_arc = np.append(sample_arc, total_len)
    sample_xyz = np.column_stack([np.interp(sample_arc, arc, xyz_mm[:, axis]) for axis in range(3)])
    sample_radii = np.interp(sample_arc, arc, radii_mm)
    return sample_xyz, sample_radii, sample_arc


def detect_branch_lesions(branch: dict, case_id: str, side: str) -> List[dict]:
    sample_xyz, sample_radii, sample_arc = sample_branch_profile(branch["xyz_mm"], branch["radii_mm"], step_mm=TARGET_SPACING_MM)
    branch["sample_xyz_mm"] = sample_xyz
    branch["sample_radii_mm"] = sample_radii
    branch["sample_arc_mm"] = sample_arc
    branch["sample_step_mm"] = TARGET_SPACING_MM

    if len(sample_radii) < 3:
        return []

    candidates = detect_stenosis_candidates(
        radii_mm=sample_radii,
        step_mm=TARGET_SPACING_MM,
        stenosis_threshold=0.30,
        min_length_mm=2.0,
        max_length_mm=20.0,
        merge_gap_mm=3.0,
        reference_window_mm=5.0,
    )

    lesions = []
    for idx, lesion in enumerate(candidates, start=1):
        start_idx = lesion["start_idx"]
        end_idx = lesion["end_idx"]
        center_idx = lesion["center_idx"]
        center_xyz = sample_xyz[center_idx]
        lesion_id = f"{case_id}_{side}_b{branch['branch_id']}_l{idx}"
        lesion_points = sample_xyz[start_idx:end_idx + 1]
        lesions.append(
            {
                "case_id": case_id,
                "side": side,
                "lesion_id": lesion_id,
                "branch_id": int(branch["branch_id"]),
                "branch_level": int(branch["branch_level"]),
                "center_xyz_mm": tuple(float(v) for v in center_xyz),
                "lesion_length_mm": float(lesion["lesion_length_mm"]),
                "min_radius_mm": float(lesion["min_radius_mm"]),
                "ref_radius_mm": float(lesion["ref_radius_mm"]),
                "stenosis_pct": float(lesion["stenosis_pct"]),
                "severity_score": float(lesion["severity_score"]),
                "start_idx": int(start_idx),
                "end_idx": int(end_idx),
                "center_idx": int(center_idx),
                "branch_sample_xyz_mm": lesion_points,
                "branch_sample_radii_mm": sample_radii[start_idx:end_idx + 1],
                "full_branch_sample_xyz_mm": sample_xyz,
                "full_branch_sample_radii_mm": sample_radii,
                "full_branch_stenosis_pct": lesion["point_stenosis_pct"],
            }
        )
    return lesions


def _rotation_matrix(azimuth_deg: float, elevation_deg: float) -> np.ndarray:
    az = math.radians(azimuth_deg)
    el = math.radians(elevation_deg)
    cos_az, sin_az = math.cos(az), math.sin(az)
    cos_el, sin_el = math.cos(el), math.sin(el)
    rz = np.array([[cos_az, -sin_az, 0.0], [sin_az, cos_az, 0.0], [0.0, 0.0, 1.0]], dtype=float)
    rx = np.array([[1.0, 0.0, 0.0], [0.0, cos_el, -sin_el], [0.0, sin_el, cos_el]], dtype=float)
    return rx @ rz


def _project_points(points_mm: np.ndarray, azimuth_deg: float, elevation_deg: float, center_mm: np.ndarray) -> np.ndarray:
    rotated = (points_mm - center_mm) @ _rotation_matrix(azimuth_deg, elevation_deg).T
    return rotated[:, :2]


def _polyline_length(points_2d: np.ndarray) -> float:
    if len(points_2d) < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(points_2d, axis=0), axis=1).sum())


def _quantize_pixels(points_2d: np.ndarray, pixel_size_mm: float = TARGET_SPACING_MM) -> set:
    if len(points_2d) == 0:
        return set()
    quantized = np.rint(points_2d / float(pixel_size_mm)).astype(int)
    return {tuple(int(v) for v in row) for row in quantized}


def score_view_for_lesion(lesion: dict, all_branch_samples: List[dict], azimuth_deg: float, elevation_deg: float) -> dict:
    lesion_points = lesion["branch_sample_xyz_mm"]
    center_mm = np.asarray(lesion["center_xyz_mm"], dtype=float)
    lesion_proj = _project_points(lesion_points, azimuth_deg, elevation_deg, center_mm)
    real_length = max(float(lesion["lesion_length_mm"]), 1e-6)
    proj_length = _polyline_length(lesion_proj)
    foreshortening = min(1.0, proj_length / real_length)

    other_points = []
    lesion_id = lesion["lesion_id"]
    for sample in all_branch_samples:
        if sample["lesion_id"] == lesion_id:
            continue
        other_points.append(sample["points_mm"])
    other_points = np.vstack(other_points) if other_points else np.empty((0, 3), dtype=float)
    other_proj = _project_points(other_points, azimuth_deg, elevation_deg, center_mm) if len(other_points) else np.empty((0, 2), dtype=float)

    lesion_pixels = _quantize_pixels(lesion_proj)
    other_pixels = _quantize_pixels(other_proj)
    overlap = 1.0
    if lesion_pixels:
        overlap = len(lesion_pixels - other_pixels) / float(len(lesion_pixels))

    total = 0.6 * foreshortening + 0.4 * overlap
    return {
        "case_id": lesion["case_id"],
        "side": lesion["side"],
        "target": "local",
        "lesion_id_or_all": lesion["lesion_id"],
        "azimuth_deg": int(azimuth_deg),
        "elevation_deg": int(elevation_deg),
        "total_score": float(total),
        "foreshortening_score": float(foreshortening),
        "overlap_score": float(overlap),
        "covered_lesions": "",
    }


def evaluate_local_views(lesion: dict, all_branch_samples: List[dict]):
    all_scores = []
    for azimuth_deg in AZIMUTH_GRID:
        for elevation_deg in ELEVATION_GRID:
            all_scores.append(score_view_for_lesion(lesion, all_branch_samples, azimuth_deg, elevation_deg))
    top_scores = select_top_local_views(all_scores, top_k=LOCAL_TOP_K, min_angle_deg=15.0)
    for rank, row in enumerate(top_scores, start=1):
        row["rank"] = rank
    return all_scores, top_scores


def build_global_view_rows(
    score_map: Dict[str, List[dict]],
    case_id: str,
    side: str,
    max_views_list: Iterable[int] = (2, 3),
) -> List[dict]:
    lesion_ids = list(score_map.keys())
    if not lesion_ids:
        return []

    best_scores = {lesion_id: max(row["total_score"] for row in rows) for lesion_id, rows in score_map.items()}
    view_buckets = {}
    for lesion_id, rows in score_map.items():
        threshold = 0.7 * best_scores[lesion_id]
        for row in rows:
            if row["total_score"] + 1e-9 < threshold:
                continue
            key = (row["azimuth_deg"], row["elevation_deg"])
            bucket = view_buckets.setdefault(
                key,
                {
                    "view_id": f"az{row['azimuth_deg']}_el{row['elevation_deg']}",
                    "azimuth_deg": row["azimuth_deg"],
                    "elevation_deg": row["elevation_deg"],
                    "covered_lesions": set(),
                    "scores": [],
                },
            )
            bucket["covered_lesions"].add(lesion_id)
            bucket["scores"].append(row["total_score"])

    candidates = []
    for bucket in view_buckets.values():
        candidates.append(
            {
                "view_id": bucket["view_id"],
                "azimuth_deg": bucket["azimuth_deg"],
                "elevation_deg": bucket["elevation_deg"],
                "covered_lesions": set(bucket["covered_lesions"]),
                "mean_score": float(np.mean(bucket["scores"])) if bucket["scores"] else 0.0,
            }
        )

    rows = []
    for max_views in max_views_list:
        selected = greedy_cover_views(lesion_ids, candidates, max_views=max_views)
        for rank, item in enumerate(selected, start=1):
            rows.append(
                {
                    "case_id": case_id,
                    "side": side,
                    "target": f"global_{max_views}view",
                    "lesion_id_or_all": "all",
                    "azimuth_deg": int(item["azimuth_deg"]),
                    "elevation_deg": int(item["elevation_deg"]),
                    "total_score": float(item.get("mean_score", 0.0)),
                    "foreshortening_score": np.nan,
                    "overlap_score": np.nan,
                    "covered_lesions": ";".join(sorted(item["covered_lesions"])),
                    "rank": rank,
                }
            )
    return rows


def analyze_coronary_mask(mask: np.ndarray, spacing: Tuple[float, float, float], case_id: str, side: str) -> dict:
    processed_mask, processed_spacing, prep_stats = preprocess_mask(mask, spacing, target_spacing=TARGET_SPACING_MM)
    radius_map = exact_edt(processed_mask, processed_spacing)
    centerline_graph, centerline_meta = extract_centerline_graph(processed_mask, processed_spacing)
    branches, root_node, aux = decompose_branches(centerline_graph, radius_map, processed_spacing)

    lesions = []
    for branch in branches:
        lesions.extend(detect_branch_lesions(branch, case_id, side))
    lesions = sorted(lesions, key=lambda item: item["severity_score"], reverse=True)
    primary_lesions = lesions[:3]

    all_branch_samples = []
    for branch in branches:
        sample_xyz = branch.get("sample_xyz_mm")
        if sample_xyz is None:
            sample_xyz, sample_radii, sample_arc = sample_branch_profile(branch["xyz_mm"], branch["radii_mm"], step_mm=TARGET_SPACING_MM)
            branch["sample_xyz_mm"] = sample_xyz
            branch["sample_radii_mm"] = sample_radii
            branch["sample_arc_mm"] = sample_arc
        all_branch_samples.append(
            {
                "branch_id": branch["branch_id"],
                "lesion_id": "",
                "points_mm": branch["sample_xyz_mm"],
            }
        )

    score_map = {}
    local_rows = []
    for lesion in primary_lesions:
        lesion_entry = {
            "branch_id": lesion["branch_id"],
            "lesion_id": lesion["lesion_id"],
            "points_mm": lesion["branch_sample_xyz_mm"],
        }
        view_samples = [item for item in all_branch_samples if item["branch_id"] != lesion["branch_id"]]
        view_samples.append(lesion_entry)
        all_scores, top_scores = evaluate_local_views(lesion, view_samples)
        score_map[lesion["lesion_id"]] = all_scores
        local_rows.extend(top_scores)

    global_rows = build_global_view_rows(score_map, case_id=case_id, side=side)

    lesion_rows = []
    for lesion in primary_lesions:
        lesion_rows.append(
            {
                "case_id": lesion["case_id"],
                "side": lesion["side"],
                "lesion_id": lesion["lesion_id"],
                "center_xyz_mm": ",".join(f"{v:.3f}" for v in lesion["center_xyz_mm"]),
                "branch_id": lesion["branch_id"],
                "branch_level": lesion["branch_level"],
                "lesion_length_mm": lesion["lesion_length_mm"],
                "min_radius_mm": lesion["min_radius_mm"],
                "ref_radius_mm": lesion["ref_radius_mm"],
                "stenosis_pct": lesion["stenosis_pct"],
                "severity_score": lesion["severity_score"],
            }
        )

    basic_stats = {
        "case_id": case_id,
        "side": side,
        "raw_shape": "x".join(str(v) for v in prep_stats["raw_shape"]),
        "raw_spacing_mm": ",".join(f"{v:.3f}" for v in prep_stats["raw_spacing_mm"]),
        "processed_shape": "x".join(str(v) for v in prep_stats["processed_shape"]),
        "processed_spacing_mm": ",".join(f"{v:.3f}" for v in prep_stats["processed_spacing_mm"]),
        "raw_volume_mm3": prep_stats["raw_volume_mm3"],
        "processed_volume_mm3": prep_stats["processed_volume_mm3"],
        "volume_delta_pct": prep_stats["volume_delta_pct"],
        "centerline_nodes": int(centerline_graph.number_of_nodes()),
        "centerline_edges": int(centerline_graph.number_of_edges()),
        "branch_count": int(len(branches)),
        "root_node": str(root_node),
        "lesion_count": int(len(primary_lesions)),
    }

    return {
        "case_id": case_id,
        "side": side,
        "raw_mask": mask,
        "raw_spacing": spacing,
        "processed_mask": processed_mask,
        "processed_spacing": processed_spacing,
        "radius_map": radius_map,
        "centerline_graph": centerline_graph,
        "centerline_meta": centerline_meta,
        "branches": branches,
        "root_node": root_node,
        "branch_aux": aux,
        "lesions": primary_lesions,
        "all_lesions": lesions,
        "score_map": score_map,
        "basic_stats": basic_stats,
        "lesion_rows": lesion_rows,
        "local_view_rows": local_rows,
        "global_view_rows": global_rows,
    }


def analyze_mask_path(path: Path) -> dict:
    case_id = path.parts[-3]
    side = "left" if "left" in path.parts[-2] else "right"
    mask, spacing, _ = read_nifti_mask(path)
    return analyze_coronary_mask(mask=mask, spacing=spacing, case_id=case_id, side=side)


def analyze_dataset(root: str | Path) -> dict:
    root = Path(root)
    paths = sorted(root.glob("case*/**/mask.nii.gz"))
    results = [analyze_mask_path(path) for path in paths]

    case_global_rows = []
    sensitivity_rows = []

    for result in results:
        for branch in result["branches"]:
            if len(branch.get("sample_radii_mm", [])) < 3:
                continue
            for threshold in (0.25, 0.30, 0.35):
                lesions = detect_stenosis_candidates(
                    branch["sample_radii_mm"],
                    step_mm=TARGET_SPACING_MM,
                    stenosis_threshold=threshold,
                    min_length_mm=2.0,
                    max_length_mm=20.0,
                    merge_gap_mm=3.0,
                    reference_window_mm=5.0,
                )
                top = max(lesions, key=lambda item: item["severity_score"]) if lesions else None
                sensitivity_rows.append(
                    {
                        "case_id": result["case_id"],
                        "side": result["side"],
                        "branch_id": branch["branch_id"],
                        "threshold": threshold,
                        "top_center_idx": top["center_idx"] if top else np.nan,
                        "top_stenosis_pct": top["stenosis_pct"] if top else np.nan,
                        "lesion_count": len(lesions),
                    }
                )

    cases = sorted({result["case_id"] for result in results})
    for case_id in cases:
        case_results = [result for result in results if result["case_id"] == case_id]
        merged_score_map = {}
        for result in case_results:
            merged_score_map.update(result["score_map"])
        case_global_rows.extend(build_global_view_rows(merged_score_map, case_id=case_id, side="case"))

    return {
        "results": results,
        "case_global_rows": case_global_rows,
        "sensitivity_rows": sensitivity_rows,
    }


def _plot_workflow_diagram(out_path: Path):
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.axis("off")
    xs = [0.07, 0.31, 0.55, 0.79]
    labels = ["Mask preprocess", "Centerline + radii", "Stenosis scoring", "View optimization"]
    for x, label in zip(xs, labels):
        ax.add_patch(plt.Rectangle((x - 0.1, 0.35), 0.2, 0.28, edgecolor="#2B3A67", facecolor="#E8EEF9", lw=2))
        ax.text(x, 0.49, label, ha="center", va="center", fontsize=11)
    for left, right in zip(xs[:-1], xs[1:]):
        ax.annotate("", xy=(right - 0.11, 0.49), xytext=(left + 0.11, 0.49), arrowprops=dict(arrowstyle="->", lw=2))
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _plot_mip_overview(results: List[dict], out_path: Path):
    fig, axes = plt.subplots(len(results), 3, figsize=(9, 3 * len(results)))
    if len(results) == 1:
        axes = np.asarray([axes])
    for row, result in enumerate(results):
        mask = result["processed_mask"]
        views = [
            ("XY", mask.max(axis=2).T),
            ("XZ", mask.max(axis=1).T),
            ("YZ", mask.max(axis=0).T),
        ]
        for col, (title, image) in enumerate(views):
            ax = axes[row, col]
            ax.imshow(image, cmap="gray")
            ax.set_title(f"{result['case_id']} {result['side']} {title}")
            ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _plot_centerline_overview(results: List[dict], out_path: Path):
    fig, axes = plt.subplots(len(results), 3, figsize=(9, 3 * len(results)))
    if len(results) == 1:
        axes = np.asarray([axes])
    for row, result in enumerate(results):
        mask = result["processed_mask"]
        branch_points = np.vstack([branch["sample_xyz_mm"] for branch in result["branches"] if len(branch["sample_xyz_mm"])]) if result["branches"] else np.empty((0, 3))
        processed_spacing = np.asarray(result["processed_spacing"], dtype=float)
        branch_ijk = branch_points / processed_spacing if len(branch_points) else np.empty((0, 3))
        view_data = [
            ("XY", mask.max(axis=2).T, branch_ijk[:, [0, 1]] if len(branch_ijk) else np.empty((0, 2))),
            ("XZ", mask.max(axis=1).T, branch_ijk[:, [0, 2]] if len(branch_ijk) else np.empty((0, 2))),
            ("YZ", mask.max(axis=0).T, branch_ijk[:, [1, 2]] if len(branch_ijk) else np.empty((0, 2))),
        ]
        for col, (title, image, pts) in enumerate(view_data):
            ax = axes[row, col]
            ax.imshow(image, cmap="gray")
            if len(pts):
                ax.scatter(pts[:, 0], pts[:, 1], s=1, c="#D62828")
            ax.set_title(f"{result['case_id']} {result['side']} {title}")
            ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _plot_radius_curves(results: List[dict], out_path: Path):
    fig, axes = plt.subplots(len(results), 1, figsize=(10, 2.6 * len(results)))
    if len(results) == 1:
        axes = [axes]
    for ax, result in zip(axes, results):
        top = result["lesions"][0] if result["lesions"] else None
        if top is None:
            ax.text(0.5, 0.5, f"{result['case_id']} {result['side']}: no major lesion", ha="center", va="center")
            ax.set_axis_off()
            continue
        branch = next(branch for branch in result["branches"] if branch["branch_id"] == top["branch_id"])
        x = branch["sample_arc_mm"]
        y = branch["sample_radii_mm"]
        sten = top["full_branch_stenosis_pct"]
        ax.plot(x, y, color="#1D3557", lw=1.8, label="radius (mm)")
        ax.axvspan(x[top["start_idx"]], x[top["end_idx"]], color="#E63946", alpha=0.18, label="lesion")
        ax2 = ax.twinx()
        ax2.plot(x[: len(sten)], sten, color="#E63946", ls="--", lw=1.2, label="stenosis")
        ax.set_title(f"{result['case_id']} {result['side']} | {top['lesion_id']} | stenosis={top['stenosis_pct']:.2f}")
        ax.set_xlabel("Arc length (mm)")
        ax.set_ylabel("Radius (mm)")
        ax2.set_ylabel("Stenosis")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _plot_local_view_comparison(results: List[dict], out_path: Path):
    fig, axes = plt.subplots(len(results), 3, figsize=(9, 3 * len(results)))
    if len(results) == 1:
        axes = np.asarray([axes])
    for row, result in enumerate(results):
        top = result["lesions"][0] if result["lesions"] else None
        local_rows = [row_item for row_item in result["local_view_rows"] if row_item.get("lesion_id_or_all") == (top["lesion_id"] if top else None)]
        if top is None or not local_rows:
            for col in range(3):
                ax = axes[row, col]
                ax.text(0.5, 0.5, f"{result['case_id']} {result['side']}\nno local lesion view", ha="center", va="center")
                ax.axis("off")
            continue
        branch = next(branch for branch in result["branches"] if branch["branch_id"] == top["branch_id"])
        points = branch["sample_xyz_mm"]
        center = np.asarray(top["center_xyz_mm"], dtype=float)
        for col, row_item in enumerate(local_rows[:3]):
            ax = axes[row, col]
            proj = _project_points(points, row_item["azimuth_deg"], row_item["elevation_deg"], center)
            ax.plot(proj[:, 0], proj[:, 1], color="#457B9D", lw=1.2)
            lesion_proj = _project_points(top["branch_sample_xyz_mm"], row_item["azimuth_deg"], row_item["elevation_deg"], center)
            ax.plot(lesion_proj[:, 0], lesion_proj[:, 1], color="#E63946", lw=2.2)
            ax.set_title(
                f"{result['case_id']} {result['side']}\naz={row_item['azimuth_deg']} el={row_item['elevation_deg']}\nscore={row_item['total_score']:.2f}"
            )
            ax.set_aspect("equal", adjustable="box")
            ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _plot_global_heatmap(dataset: dict, out_path: Path):
    results = dataset["results"]
    case_ids = sorted({result["case_id"] for result in results})
    fig, axes = plt.subplots(len(case_ids), 1, figsize=(10, 3.2 * len(case_ids)))
    if len(case_ids) == 1:
        axes = [axes]
    for ax, case_id in zip(axes, case_ids):
        case_results = [result for result in results if result["case_id"] == case_id]
        mean_scores = np.zeros((len(ELEVATION_GRID), len(AZIMUTH_GRID)), dtype=float)
        counts = np.zeros_like(mean_scores)
        for result in case_results:
            for rows in result["score_map"].values():
                for row in rows:
                    el_idx = ELEVATION_GRID.index(row["elevation_deg"])
                    az_idx = AZIMUTH_GRID.index(row["azimuth_deg"])
                    mean_scores[el_idx, az_idx] += row["total_score"]
                    counts[el_idx, az_idx] += 1
        with np.errstate(invalid="ignore"):
            heat = np.divide(mean_scores, counts, out=np.zeros_like(mean_scores), where=counts > 0)
        im = ax.imshow(heat, aspect="auto", cmap="viridis", origin="lower")
        ax.set_title(f"{case_id} mean local-view score")
        ax.set_xlabel("Azimuth index (5° step)")
        ax.set_ylabel("Elevation index (5° step)")
        fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def write_analysis_outputs(dataset: dict, out_dir: str | Path):
    out_dir = Path(out_dir)
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    basic_rows = [result["basic_stats"] for result in dataset["results"]]
    lesion_rows = [row for result in dataset["results"] for row in result["lesion_rows"]]
    local_rows = [row for result in dataset["results"] for row in result["local_view_rows"]]
    global_rows = [row for result in dataset["results"] for row in result["global_view_rows"]] + list(dataset.get("case_global_rows", []))
    sensitivity_rows = list(dataset.get("sensitivity_rows", []))

    pd.DataFrame(basic_rows).to_csv(tables_dir / "basic_structure_stats.csv", index=False)
    pd.DataFrame(lesion_rows).to_csv(tables_dir / "lesion_candidates.csv", index=False)
    pd.DataFrame(local_rows).to_csv(tables_dir / "local_best_views.csv", index=False)
    pd.DataFrame(global_rows).to_csv(tables_dir / "global_view_plans.csv", index=False)
    pd.DataFrame(sensitivity_rows).to_csv(tables_dir / "sensitivity_analysis.csv", index=False)

    _plot_workflow_diagram(figures_dir / "figure1_workflow.png")
    if dataset["results"] and all("processed_mask" in result for result in dataset["results"]):
        _plot_mip_overview(dataset["results"], figures_dir / "figure2_mip_overview.png")
        _plot_centerline_overview(dataset["results"], figures_dir / "figure3_centerline_overview.png")
        _plot_radius_curves(dataset["results"], figures_dir / "figure4_radius_curves.png")
        _plot_local_view_comparison(dataset["results"], figures_dir / "figure5_local_view_comparison.png")
        _plot_global_heatmap(dataset, figures_dir / "figure6_global_heatmap.png")

    summary_path = out_dir / "summary.md"
    summary_lines = [
        "# Coronary Stenosis Analysis Summary",
        "",
        f"- analyzed_masks: {len(dataset['results'])}",
        f"- table_dir: {tables_dir}",
        f"- figure_dir: {figures_dir}",
    ]
    summary_path.write_text("\n".join(summary_lines), encoding="utf-8")


def run_full_analysis(root: str | Path, out_dir: str | Path):
    dataset = analyze_dataset(root)
    write_analysis_outputs(dataset, out_dir)
    return dataset
