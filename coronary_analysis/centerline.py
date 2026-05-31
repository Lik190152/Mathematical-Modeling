from __future__ import annotations

from collections import deque
from itertools import product
from typing import Dict, Iterable, List, Tuple

import networkx as nx
import numpy as np
from skimage.morphology import skeletonize_3d

from .geometry import voxel_to_world


_NEIGHBOR_OFFSETS = [
    np.array(offset, dtype=int)
    for offset in product((-1, 0, 1), repeat=3)
    if offset != (0, 0, 0)
]


def skeletonize_volume(mask: np.ndarray) -> np.ndarray:
    return skeletonize_3d(mask.astype(np.uint8)) > 0


def extract_centerline_graph(mask: np.ndarray, spacing: Tuple[float, float, float]) -> tuple[nx.Graph, Dict[str, object]]:
    skeleton = skeletonize_volume(mask)
    coords = np.argwhere(skeleton)
    graph = nx.Graph()
    if coords.size == 0:
        return graph, {"node_count": 0, "edge_count": 0, "skeleton": skeleton}

    coord_set = {tuple(int(v) for v in coord) for coord in coords}
    spacing_arr = np.asarray(spacing, dtype=float)

    for coord in coord_set:
        graph.add_node(coord, coord_ijk=np.asarray(coord, dtype=int), xyz_mm=np.asarray(coord, dtype=float) * spacing_arr)

    for coord in coord_set:
        base = np.asarray(coord, dtype=int)
        for offset in _NEIGHBOR_OFFSETS:
            neigh = tuple(int(v) for v in (base + offset))
            if neigh not in coord_set:
                continue
            if neigh <= coord:
                continue
            weight = float(np.linalg.norm(offset * spacing_arr))
            graph.add_edge(coord, neigh, weight=weight)

    meta = {
        "node_count": graph.number_of_nodes(),
        "edge_count": graph.number_of_edges(),
        "skeleton": skeleton,
    }
    return graph, meta


def choose_root_node(graph: nx.Graph, radius_map: np.ndarray) -> tuple[int, int, int]:
    if graph.number_of_nodes() == 0:
        raise ValueError("Cannot choose root from empty graph")
    endpoints = [node for node, degree in graph.degree() if degree == 1]
    candidates = endpoints or list(graph.nodes)
    return max(candidates, key=lambda node: float(radius_map[node]))


def decompose_branches(
    graph: nx.Graph,
    radius_map: np.ndarray,
    spacing: Tuple[float, float, float],
) -> tuple[list[dict], tuple[int, int, int], dict]:
    if graph.number_of_nodes() == 0:
        return [], None, {}

    root = choose_root_node(graph, radius_map)
    branch_nodes = {node for node, degree in graph.degree() if degree != 2}
    if not branch_nodes:
        # Straight single branch: use endpoints as artificial branch nodes.
        endpoints = [node for node, degree in graph.degree() if degree == 1]
        branch_nodes = set(endpoints[:2] or [next(iter(graph.nodes))])

    visited_edges = set()
    branches = []

    for start in list(branch_nodes):
        for nbr in graph.neighbors(start):
            edge_key = frozenset((start, nbr))
            if edge_key in visited_edges:
                continue

            path = [start, nbr]
            visited_edges.add(edge_key)
            prev = start
            current = nbr
            while current not in branch_nodes:
                next_nodes = [node for node in graph.neighbors(current) if node != prev]
                if not next_nodes:
                    break
                nxt = next_nodes[0]
                visited_edges.add(frozenset((current, nxt)))
                path.append(nxt)
                prev, current = current, nxt

            coords = np.asarray(path, dtype=int)
            xyz_mm = voxel_to_world(coords, spacing)
            step_mm = np.linalg.norm(np.diff(xyz_mm, axis=0), axis=1) if len(path) > 1 else np.array([], dtype=float)
            length_mm = float(step_mm.sum())
            branch = {
                "branch_id": len(branches),
                "start_node": path[0],
                "end_node": path[-1],
                "path_nodes": path,
                "coords_ijk": coords,
                "xyz_mm": xyz_mm,
                "radii_mm": np.asarray([float(radius_map[node]) for node in path], dtype=float),
                "length_mm": length_mm,
            }
            branches.append(branch)

    branch_graph = nx.Graph()
    for branch in branches:
        branch_graph.add_node(branch["branch_id"])

    node_to_branch_ids: Dict[tuple[int, int, int], list[int]] = {}
    for branch in branches:
        for node in (branch["start_node"], branch["end_node"]):
            node_to_branch_ids.setdefault(node, []).append(branch["branch_id"])
    for branch_ids in node_to_branch_ids.values():
        for i in range(len(branch_ids)):
            for j in range(i + 1, len(branch_ids)):
                branch_graph.add_edge(branch_ids[i], branch_ids[j])

    root_branch = max(
        range(len(branches)),
        key=lambda idx: max(float(radius_map[branches[idx]["start_node"]]), float(radius_map[branches[idx]["end_node"]])),
    )
    levels = nx.single_source_shortest_path_length(branch_graph, root_branch) if branch_graph.number_of_nodes() else {}
    for branch in branches:
        branch["branch_level"] = int(levels.get(branch["branch_id"], 0))

    aux = {
        "root_node": root,
        "root_branch_id": root_branch,
        "branch_graph": branch_graph,
    }
    return branches, root, aux

