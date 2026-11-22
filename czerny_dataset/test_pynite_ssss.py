"""Regression checks for the executed SSSS Rect benchmark."""

import io
import unittest
from contextlib import redirect_stdout
from math import fsum
from types import SimpleNamespace

import numpy as np

from scripts.run_pynite_ssss import center_displacement, reaction_balance, response, solve_mesh


class TestPyNiteSSSS(unittest.TestCase):
    def test_signed_reaction_residual(self):
        balanced = reaction_balance(-100, [57, 48, -5])
        self.assertEqual(balanced["support_fz_kN"], 100)
        self.assertEqual(balanced["residual_fz_kN"], 0)
        self.assertEqual(balanced["residual_percent"], 0)
        self.assertEqual(reaction_balance(-100, [102])["residual_percent"], 2)
        self.assertEqual(reaction_balance(-100, [-100])["residual_percent"], -200)

    def test_moment_denominator_and_deflection_multiplier(self):
        self.assertEqual(response(8, 160, denominator=True)["coefficient"], 20)
        self.assertEqual(response(0.002, 0.04, denominator=False)["coefficient"], 0.05)

    def test_center_displacement_uses_the_rect_bending_polynomial(self):
        # A combination of basis terms, with rotations from the local slopes.
        def w(x, y):
            return x**3 * y + 2 * x * y**3 + x**2 + y**2

        dx, dy = 0.8, 1.2
        nodal_dofs = np.zeros(24)
        for i, (x, y) in enumerate(((0, 0), (dx, 0), (dx, dy), (0, dy))):
            rx = x**3 + 6 * x * y**2 + 2 * y
            ry = -(3 * x**2 * y + 2 * y**3 + 2 * x)
            nodal_dofs[6 * i + 2:6 * i + 5] = w(x, y), rx, ry
        plate = SimpleNamespace(d=lambda _: nodal_dofs, width=lambda: dx, height=lambda: dy)
        self.assertAlmostEqual(center_displacement(plate), w(dx / 2, dy / 2))

    def test_rectangular_solver_reports_off_center_sagging_peak_and_signed_supports(self):
        with redirect_stdout(io.StringIO()):
            result = solve_mesh(ratio=2.0, nx=7)
        self.assertEqual(result["center_element_local_axes_in_global"], np.eye(3).tolist())
        reactions = result["reactions"]
        self.assertAlmostEqual(reactions["support_fz_kN"], 320, places=6)
        self.assertAlmostEqual(reactions["residual_percent"], 0, places=6)
        positions = [node for group in reactions["support_nodes_by_location"].values() for node in group]
        self.assertEqual(len({(node["x_m"], node["y_m"]) for node in positions}), 2 * (7 + 15))
        corners = reactions["support_nodes_by_location"]["corners"]
        self.assertTrue(all(node["fz_kN"] < 0 for node in corners))
        self.assertAlmostEqual(
            fsum(reactions["open_edge_resultants_kN"].values()) + fsum(node["fz_kN"] for node in corners),
            reactions["support_fz_kN"],
        )
        center = result["center"]
        self.assertGreater(center["Mx"], center["My"])
        self.assertLess(center["w_m"], 0)
        self.assertAlmostEqual(center["Mxy"], 0, places=9)
        self.assertAlmostEqual(result["responses"]["Mx_mid"]["coefficient"], 10 * 4**2 / center["Mx"])
        self.assertAlmostEqual(
            result["responses"]["downward_deflection_mid"]["coefficient"],
            -center["w_m"] * 33e6 * 0.18**3 / (10 * 4**4),
        )
        extrema = result["positive_moment_extrema"]
        self.assertGreater(extrema["My_global"]["My"], 1.25 * center["My"])
        self.assertNotAlmostEqual(extrema["My_mid_section"]["y_m"], center["y_m"])
        self.assertEqual(extrema["My_mid_section"]["x_m"], center["x_m"])


if __name__ == "__main__":
    unittest.main()
