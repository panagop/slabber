"""Regression tests for independent SSSS FD solver + PyNite model setup helpers.

These do not certify the original Czerny transcription or run PyNite.
"""
import unittest
from finite_difference_reference import solve_ssss
from pynite_czerny_benchmark import edge_codes, nearest_odd

class TestFDReference(unittest.TestCase):
    def test_square_center_coefficients(self):
        v=solve_ssss(nx=32,ny=32)
        self.assertLess(abs(v['kx']/27.2-1),.01)
        self.assertLess(abs(v['ky']/27.2-1),.01)
        self.assertLess(abs(v['deflection_multiplier']/0.0487-1),.01)

    def test_rectangular_my_max_is_not_at_center(self):
        v=solve_ssss(b=8,nx=32,ny=64)
        self.assertGreater(v['My_max'],v['My_mid']*1.30)
        self.assertLess(abs(v['ky']/40.3-1),.01)
        self.assertLess(abs(v['kx']/10.4-1),.01)

    def test_source_edge_order(self):
        self.assertEqual(edge_codes('SSSS'),['S','S','S','S'])
        self.assertEqual(edge_codes('H17_SSS'),['S','FREE','S','S'])
        self.assertEqual(edge_codes('H22_ALL_FIXED'),['F','FREE','F','F'])

    def test_mesh_centering(self):
        self.assertEqual(nearest_odd(11),11)
        self.assertEqual(nearest_odd(12),13)

if __name__=='__main__':
    unittest.main()
