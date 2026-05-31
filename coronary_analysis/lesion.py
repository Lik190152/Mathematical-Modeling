from __future__ import annotations

import math
from typing import List

import numpy as np


def smooth_profile(values: np.ndarray, window: int = 5) -> np.ndarray:
    if len(values) == 0:
        return values.astype(float)
    if window <= 1:
        return values.astype(float)
    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(values, kernel, mode="same")


def detect_stenosis_candidates(
    radii_mm: np.ndarray,
    step_mm: float,
    stenosis_threshold: float = 0.30,
    min_length_mm: float = 2.0,
    max_length_mm: float = 20.0,
    merge_gap_mm: float = 3.0,
    reference_window_mm: float = 5.0,
) -> List[dict]:
    radii_mm = np.asarray(radii_mm, dtype=float)
    if radii_mm.size == 0:
        return []

    smooth = smooth_profile(radii_mm, window=5)
    window_pts = max(1, int(round(reference_window_mm / step_mm)))
    stenosis = np.zeros_like(smooth)
    ref_radius = np.zeros_like(smooth)

    for idx in range(len(smooth)):
        before = smooth[max(0, idx - window_pts):idx]
        after = smooth[idx + 1:min(len(smooth), idx + 1 + window_pts)]
        if before.size == 0 or after.size == 0:
            continue
        prox = float(before.mean())
        dist = float(after.mean())
        ref = 0.5 * (prox + dist)
        ref_radius[idx] = ref
        if ref > 1e-6:
            stenosis[idx] = max(0.0, 1.0 - (2.0 * smooth[idx] / (prox + dist)))

    active = stenosis >= stenosis_threshold
    intervals = []
    start = None
    for idx, flag in enumerate(active):
        if flag and start is None:
            start = idx
        elif not flag and start is not None:
            intervals.append([start, idx - 1])
            start = None
    if start is not None:
        intervals.append([start, len(active) - 1])

    merged = []
    for interval in intervals:
        if not merged:
            merged.append(interval)
            continue
        gap_mm = (interval[0] - merged[-1][1] - 1) * step_mm
        if gap_mm < merge_gap_mm:
            merged[-1][1] = interval[1]
        else:
            merged.append(interval)

    lesions = []
    for start, end in merged:
        length_mm = (end - start + 1) * step_mm
        if length_mm < min_length_mm or length_mm > max_length_mm:
            continue
        local_idx = start + int(np.argmin(smooth[start:end + 1]))
        min_radius = float(smooth[local_idx])
        ref = float(ref_radius[local_idx]) if ref_radius[local_idx] > 0 else float(np.mean(smooth[start:end + 1]))
        stenosis_pct = float(stenosis[local_idx])
        severity = stenosis_pct * math.log1p(length_mm)
        lesions.append(
            {
                "start_idx": int(start),
                "end_idx": int(end),
                "center_idx": int((start + end) // 2),
                "min_idx": int(local_idx),
                "lesion_length_mm": float(length_mm),
                "min_radius_mm": min_radius,
                "ref_radius_mm": ref,
                "stenosis_pct": stenosis_pct,
                "severity_score": float(severity),
                "smoothed_radii_mm": smooth,
                "point_stenosis_pct": stenosis,
            }
        )

    return lesions

