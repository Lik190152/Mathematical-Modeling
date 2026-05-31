import tempfile
from pathlib import Path
import math
import unittest

import numpy as np

from coronary_analysis.geometry import (
    exact_edt,
    keep_largest_component,
)
from coronary_analysis.centerline import extract_centerline_graph
from coronary_analysis.lesion import detect_stenosis_candidates
from coronary_analysis.pipeline import analyze_coronary_mask, write_analysis_outputs
from coronary_analysis.view_selection import (
    greedy_cover_views,
    select_top_local_views,
)


class GeometryTests(unittest.TestCase):
    def test_keep_largest_component_discards_smaller_fragments(self):
        mask = np.zeros((6, 6, 6), dtype=bool)
        mask[1:3, 1:3, 1:3] = True
        mask[4, 4, 4] = True
        kept = keep_largest_component(mask)
        self.assertEqual(int(kept.sum()), 8)
        self.assertTrue(kept[1, 1, 1])
        self.assertFalse(kept[4, 4, 4])

    def test_exact_edt_matches_simple_cube_distances(self):
        mask = np.zeros((5, 5, 5), dtype=bool)
        mask[1:4, 1:4, 1:4] = True
        dist = exact_edt(mask, spacing=(1.0, 1.0, 1.0))
        self.assertAlmostEqual(float(dist[2, 2, 2]), 2.0, places=6)
        self.assertAlmostEqual(float(dist[1, 2, 2]), 1.0, places=6)
        self.assertEqual(float(dist[0, 0, 0]), 0.0)

    def test_extract_centerline_graph_recovers_t_shape_topology(self):
        mask = np.zeros((7, 7, 7), dtype=bool)
        mask[3, 3, 1:6] = True
        mask[3:6, 3, 3] = True
        graph, meta = extract_centerline_graph(mask, spacing=(1.0, 1.0, 1.0))
        degrees = sorted(dict(graph.degree()).values())
        self.assertEqual(degrees.count(1), 3)
        self.assertEqual(degrees.count(3), 1)
        self.assertEqual(meta["node_count"], graph.number_of_nodes())


class LesionTests(unittest.TestCase):
    def test_detect_stenosis_candidates_finds_single_dip(self):
        step_mm = 0.4
        radii = np.full(80, 2.0, dtype=float)
        radii[25:40] = np.array([1.9, 1.7, 1.5, 1.2, 1.0, 0.9, 0.9, 0.9, 0.9, 0.9, 1.0, 1.2, 1.5, 1.7, 1.9])
        lesions = detect_stenosis_candidates(
            radii_mm=radii,
            step_mm=step_mm,
            stenosis_threshold=0.30,
            min_length_mm=2.0,
            max_length_mm=20.0,
            merge_gap_mm=3.0,
            reference_window_mm=5.0,
        )
        self.assertEqual(len(lesions), 1)
        lesion = lesions[0]
        self.assertGreaterEqual(lesion["stenosis_pct"], 0.40)
        self.assertGreaterEqual(lesion["lesion_length_mm"], 2.0)
        self.assertLessEqual(lesion["lesion_length_mm"], 20.0)


class ViewSelectionTests(unittest.TestCase):
    def test_select_top_local_views_enforces_angular_separation(self):
        candidates = [
            {"azimuth_deg": 0, "elevation_deg": 0, "total_score": 0.95},
            {"azimuth_deg": 5, "elevation_deg": 0, "total_score": 0.94},
            {"azimuth_deg": 15, "elevation_deg": 0, "total_score": 0.93},
            {"azimuth_deg": 45, "elevation_deg": 0, "total_score": 0.92},
            {"azimuth_deg": 90, "elevation_deg": 0, "total_score": 0.91},
        ]
        selected = select_top_local_views(candidates, top_k=3, min_angle_deg=15.0)
        self.assertEqual(len(selected), 3)
        self.assertEqual(selected[0]["azimuth_deg"], 0)
        self.assertEqual(selected[1]["azimuth_deg"], 15)
        self.assertEqual(selected[2]["azimuth_deg"], 45)

    def test_greedy_cover_views_prefers_multi_lesion_coverage(self):
        lesions = ["L1", "L2", "L3"]
        candidates = [
            {"view_id": "A", "covered_lesions": {"L1", "L2"}, "mean_score": 0.9},
            {"view_id": "B", "covered_lesions": {"L2"}, "mean_score": 0.7},
            {"view_id": "C", "covered_lesions": {"L3"}, "mean_score": 0.8},
            {"view_id": "D", "covered_lesions": {"L1", "L3"}, "mean_score": 0.6},
        ]
        selected = greedy_cover_views(lesions, candidates, max_views=2)
        self.assertEqual([item["view_id"] for item in selected], ["A", "C"])


class PipelineTests(unittest.TestCase):
    def test_analyze_coronary_mask_runs_end_to_end_on_simple_tree(self):
        mask = np.zeros((9, 9, 9), dtype=bool)
        mask[4, 4, 1:8] = True
        mask[4:8, 4, 4] = True
        result = analyze_coronary_mask(
            mask=mask,
            spacing=(1.0, 1.0, 1.0),
            case_id="synthetic",
            side="left",
        )
        self.assertEqual(result["basic_stats"]["case_id"], "synthetic")
        self.assertEqual(result["basic_stats"]["side"], "left")
        self.assertGreater(result["basic_stats"]["centerline_nodes"], 0)
        self.assertGreaterEqual(result["basic_stats"]["branch_count"], 1)
        self.assertIsInstance(result["lesions"], list)
        self.assertIsInstance(result["local_view_rows"], list)
        self.assertIsInstance(result["global_view_rows"], list)

    def test_write_analysis_outputs_creates_core_tables(self):
        sample = {
            "results": [
                {
                    "basic_stats": {"case_id": "caseX", "side": "left", "centerline_nodes": 10, "branch_count": 2},
                    "lesion_rows": [{"case_id": "caseX", "side": "left", "lesion_id": "L1"}],
                    "local_view_rows": [],
                    "global_view_rows": [],
                }
            ],
            "case_global_rows": [],
            "sensitivity_rows": [],
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            out_dir = Path(tmpdir) / "analysis_results"
            write_analysis_outputs(sample, out_dir)
            self.assertTrue((out_dir / "tables" / "basic_structure_stats.csv").exists())
            self.assertTrue((out_dir / "tables" / "lesion_candidates.csv").exists())


if __name__ == "__main__":
    unittest.main()
