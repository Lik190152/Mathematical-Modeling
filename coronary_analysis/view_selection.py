from __future__ import annotations

import math
from typing import Iterable, List

import numpy as np


def view_direction(azimuth_deg: float, elevation_deg: float) -> np.ndarray:
    az = math.radians(azimuth_deg)
    el = math.radians(elevation_deg)
    return np.array(
        [
            math.cos(el) * math.cos(az),
            math.cos(el) * math.sin(az),
            math.sin(el),
        ],
        dtype=float,
    )


def angular_distance_deg(view_a: dict, view_b: dict) -> float:
    vec_a = view_direction(view_a["azimuth_deg"], view_a["elevation_deg"])
    vec_b = view_direction(view_b["azimuth_deg"], view_b["elevation_deg"])
    dot = float(np.clip(np.dot(vec_a, vec_b), -1.0, 1.0))
    return math.degrees(math.acos(dot))


def select_top_local_views(candidates: List[dict], top_k: int = 3, min_angle_deg: float = 15.0) -> List[dict]:
    ordered = sorted(
        candidates,
        key=lambda item: (item["total_score"], item.get("foreshortening_score", 0.0), item.get("overlap_score", 0.0)),
        reverse=True,
    )
    selected = []
    for candidate in ordered:
        if all(angular_distance_deg(candidate, prior) + 1e-9 >= min_angle_deg for prior in selected):
            selected.append(candidate)
        if len(selected) >= top_k:
            break
    return selected


def greedy_cover_views(lesion_ids: Iterable[str], candidates: List[dict], max_views: int) -> List[dict]:
    lesion_ids = list(lesion_ids)
    remaining = set(lesion_ids)
    selected = []
    unused = list(candidates)

    while unused and len(selected) < max_views:
        best = None
        best_key = None
        for candidate in unused:
            new_cover = len(candidate["covered_lesions"] & remaining)
            key = (new_cover, candidate.get("mean_score", 0.0), candidate.get("total_score", 0.0))
            if best is None or key > best_key:
                best = candidate
                best_key = key
        selected.append(best)
        remaining -= best["covered_lesions"]
        unused = [candidate for candidate in unused if candidate is not best]
        if not remaining and len(selected) >= max_views:
            break

    return selected[:max_views]
