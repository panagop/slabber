"""Regression checks for the Hahn boundary conditions and source conventions."""

import unittest

import numpy as np

from czerny_dataset.czerny import THREE_EDGE_CSV, load_table
from czerny_dataset.finite_difference_reference import solve_ssss
from scripts.kirchhoff_ritz import solve_ritz
from scripts.run_pynite_ssss import COMBO, E, LOAD, LX, THICKNESS, solve_rect_model
from scripts.run_three_edge_benchmark import hahn_reference, negative_reach, rect_run, ritz_run, uses_rear_peak


class TestThreeEdgeBenchmark(unittest.TestCase):
    def test_ritz_matches_independent_four_edge_baseline(self):
        field = solve_ritz("SSSS", LX, 6, LOAD, E, THICKNESS, 0, 16)
        fd = solve_ssss(a=LX, b=6, nx=64, ny=96)
        mx, my, _ = field.moments(LX / 2, 3)
        np.testing.assert_allclose([mx.item(), my.item(), field.partial(LX / 2, 3).item()],
                                   [fd["Mx_mid"], fd["My_mid"], fd["w_mid"]], rtol=0.002)

    def test_hahn_scaling_and_h17_worked_example(self):
        table = load_table(THREE_EDGE_CSV)
        source = hahn_reference("H17_SSS", 0.6, table)
        load = LOAD * LX**2 * 0.6
        self.assertAlmostEqual(source["mx_free_mid"]["value"], load / 9.2)
        self.assertAlmostEqual(source["reaction_x_total"]["value"], load * 0.31)
        self.assertAlmostEqual(source["deflection_ref"]["value"], load * LX**2 / (7.4 * E * THICKNESS**3))
        with self.assertRaisesRegex(ValueError, "between printed Hahn knots"):
            hahn_reference("H17_SSS", 0.55, table)

        run = ritz_run("H17_SSS", 0.6, 20, ["S", "FREE", "S", "S"])
        for quantity in ("mx_center", "my_center", "mx_free_mid", "reaction_x_total", "reaction_y_total"):
            actual = run["points"][quantity]["value"]
            self.assertAlmostEqual(actual / source[quantity]["value"], 1, delta=0.01)
        reactions = run["reactions_from_effective_shear"]
        self.assertLess(abs(reactions["residual_percent"]), 0.05)
        forces = [point["fz_kN"] for point in reactions["signed_corner_twisting_forces"]]
        self.assertTrue(forces[0] > 0 and forces[1] > 0 and forces[2] < 0 and forces[3] < 0)
        self.assertLess(abs(run["free_edge_effective_shear_mid_kN_per_m"]), 0.01)

    def test_rect_free_edge_does_not_release_side_corner_supports(self):
        m, _ = solve_rect_model("H22_ALL_FIXED", 0.5, 7)
        free_node = m.nodes["N_3_0"]
        self.assertFalse(free_node.support_DZ or free_node.support_RX or free_node.support_RY)
        self.assertLess(free_node.DZ[COMBO], 0)
        corner = m.nodes["N_0_0"]
        self.assertTrue(corner.support_DZ and corner.support_RX and corner.support_RY)
        self.assertEqual(corner.DZ[COMBO], 0)

        result = rect_run("H17_SSS", 0.5, 7)
        self.assertAlmostEqual(result["reactions"]["support_fz_kN"], 80, places=6)
        self.assertEqual(result["reactions"]["free_edge_max_abs_nodal_fz_kN"], 0)

    def test_station_selection_and_rear_sign_change(self):
        self.assertFalse(uses_rear_peak("H19_ONE_SIDE_FIXED", 0.8))
        self.assertTrue(uses_rear_peak("H19_ONE_SIDE_FIXED", 1.0))
        self.assertFalse(uses_rear_peak("H20_BOTH_SIDES_FIXED", 0.6))
        self.assertTrue(uses_rear_peak("H20_BOTH_SIDES_FIXED", 0.8))
        reach = negative_reach([(0, 0), (1, 2), (2, -2), (3, -4)], 3)
        self.assertAlmostEqual(reach["value"], 0.5)
        self.assertAlmostEqual(reach["y_m"], 1.5)


if __name__ == "__main__":
    unittest.main()
