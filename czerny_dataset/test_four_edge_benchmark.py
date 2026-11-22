"""Regression checks for the four-edge comparison and source stations."""

import unittest

import numpy as np

from czerny_dataset.czerny import load_table
from czerny_dataset.finite_difference_reference import solve_ssss
from scripts.run_four_edge_benchmark import CASES, FOUR_EDGE_RATIOS, czerny_samples, fd_mesh, rect_mesh


class TestFourEdgeBenchmark(unittest.TestCase):
    def test_source_rows_are_exact_at_all_selectable_ratios(self):
        self.assertEqual([round(ratio * 100) for ratio in FOUR_EDGE_RATIOS], list(range(100, 201, 5)))
        table = load_table()
        for case in CASES:
            for ratio in FOUR_EDGE_RATIOS:
                with self.subTest(case=case, ratio=ratio):
                    rows = czerny_samples(case, ratio, table)
                    self.assertTrue(rows)
                    self.assertTrue(all(not row["source"]["interpolated"] for row in rows.values()))

    def test_fd_keeps_center_node_when_nominal_long_count_rounds_odd(self):
        mesh = fd_mesh(1.1, 32, "FFFF")
        self.assertEqual(mesh["ny"], 36)
        midpoint = mesh["points"]["mx_mid"]
        self.assertAlmostEqual(midpoint["x_m"], 2.0)
        self.assertAlmostEqual(midpoint["y_m"], 2.2)

    def test_simply_supported_fd_matches_inherited_reference(self):
        inherited = solve_ssss(a=4, b=6, nx=32, ny=48)
        points = fd_mesh(1.5, 32, "SSSS")["points"]
        for new, original in (
            ("mx_mid", "Mx_mid"), ("my_mid", "My_mid"), ("my_global_max", "My_max"),
            ("deflection_mid", "w_mid"),
        ):
            self.assertAlmostEqual(points[new]["value"], inherited[original], places=8)

    def test_square_rotates_single_fixed_edge(self):
        fixed_on_x = fd_mesh(1.0, 32, "FSSS")
        fixed_on_y = fd_mesh(1.0, 32, "SFSS")
        x_points, y_points = fixed_on_x["points"], fixed_on_y["points"]
        self.assertAlmostEqual(x_points["mx_mid"]["value"], y_points["my_mid"]["value"])
        self.assertAlmostEqual(x_points["mx_fixed_mid"]["value"], y_points["my_fixed_mid"]["value"])
        self.assertAlmostEqual(x_points["deflection_mid"]["value"], y_points["deflection_mid"]["value"])

    def test_adjacent_clamp_table_uses_edge_minima(self):
        source = czerny_samples("FF_adj", 1.5, load_table())
        self.assertNotIn("mx_fixed_mid", source)
        self.assertEqual(source["mx_fixed_min"]["source"]["quantity"], "mx_fixed_mid")
        self.assertEqual(source["mx_fixed_min"]["source"]["printed_symbol_from_scan"], "m_xermin")
        self.assertEqual(source["my_fixed_min"]["source"]["printed_symbol_from_scan"], "m_yermin")
        self.assertEqual(source["mx_fixed_min"]["source"]["source_page"], 286)
        self.assertAlmostEqual(source["my_fixed_min"]["value"], -160 / 12.4)
        field = fd_mesh(1.5, 32, "FFSS")["fixed_edge_moments"]
        for row in field.values():
            self.assertLess(row["minimum"]["value"], row["midpoint"]["value"])

    def test_square_fixed_rect_respects_symmetry_and_hogging(self):
        mesh = rect_mesh("FFFF", 1.0, 7, "FFFF")
        fixed = mesh["fixed_edge_moments"]
        support_moments = [row["midpoint"]["value"] for row in fixed.values()]
        np.testing.assert_allclose(support_moments, support_moments[0], atol=1e-9)
        self.assertLess(support_moments[0], 0)
        self.assertGreater(mesh["points"]["mx_mid"]["value"], 0)
        self.assertAlmostEqual(mesh["reactions"]["support_fz_kN"], 160, places=6)
        self.assertEqual(mesh["simply_supported_corner_twisting"], [])


if __name__ == "__main__":
    unittest.main()
