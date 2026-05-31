from __future__ import annotations

import gzip
import struct
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
from scipy import ndimage


def read_nifti_mask(path: str | Path) -> Tuple[np.ndarray, Tuple[float, float, float], Dict[str, float]]:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as f:
        hdr = f.read(348)
        sizeof_hdr = struct.unpack("<I", hdr[0:4])[0]
        if sizeof_hdr != 348:
            raise ValueError(f"{path}: unexpected header size {sizeof_hdr}")

        dim = struct.unpack("<8h", hdr[40:56])
        datatype = struct.unpack("<h", hdr[70:72])[0]
        bitpix = struct.unpack("<h", hdr[72:74])[0]
        pixdim = struct.unpack("<8f", hdr[76:108])
        vox_offset = int(struct.unpack("<f", hdr[108:112])[0])
        qform_code = struct.unpack("<h", hdr[252:254])[0]
        sform_code = struct.unpack("<h", hdr[254:256])[0]

        shape = tuple(dim[1:1 + dim[0]])
        spacing = tuple(float(v) for v in pixdim[1:1 + dim[0]])

        if len(shape) != 3:
            raise ValueError(f"{path}: expected 3D NIfTI, got shape {shape}")
        if datatype != 2 or bitpix != 8:
            raise ValueError(f"{path}: expected uint8 mask, got datatype={datatype}, bitpix={bitpix}")

        f.seek(vox_offset)
        raw = f.read(int(np.prod(shape)))

    mask = np.frombuffer(raw, dtype=np.uint8).reshape(shape, order="F") > 0
    meta = {
        "qform_code": float(qform_code),
        "sform_code": float(sform_code),
        "vox_offset": float(vox_offset),
    }
    return mask, spacing, meta


def keep_largest_component(mask: np.ndarray, connectivity: int = 3) -> np.ndarray:
    structure = ndimage.generate_binary_structure(mask.ndim, connectivity)
    labeled, num = ndimage.label(mask, structure=structure)
    if num <= 1:
        return mask.astype(bool, copy=True)
    counts = np.bincount(labeled.ravel())
    counts[0] = 0
    keep_label = counts.argmax()
    return labeled == keep_label


def exact_edt(mask: np.ndarray, spacing: Tuple[float, float, float]) -> np.ndarray:
    return ndimage.distance_transform_edt(mask, sampling=spacing)


def resample_binary_mask(
    mask: np.ndarray,
    spacing: Tuple[float, float, float],
    target_spacing: float = 0.4,
) -> tuple[np.ndarray, Tuple[float, float, float], np.ndarray]:
    scale = np.asarray(spacing, dtype=float) / float(target_spacing)
    resampled = ndimage.zoom(mask.astype(np.uint8), zoom=scale, order=0) > 0
    target = (float(target_spacing), float(target_spacing), float(target_spacing))
    return resampled, target, scale


def mask_volume_mm3(mask: np.ndarray, spacing: Tuple[float, float, float]) -> float:
    return float(mask.sum()) * float(np.prod(spacing))


def voxel_to_world(coords_ijk: np.ndarray, spacing: Tuple[float, float, float]) -> np.ndarray:
    spacing_arr = np.asarray(spacing, dtype=float)
    return coords_ijk.astype(float) * spacing_arr

