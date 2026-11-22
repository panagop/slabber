"""Regression checks for the 1996, printed pp. 281-289, 0.05 grid.

These checks enforce structure and a selection of printed anchors; they do NOT
replace independent visual double-key of all 1,701 transcribed numbers.
"""
import math
import unittest
from czerny import calculate, load_table, lookup


class CzernyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = load_table()

    def test_all_nine_cases_have_full_grid(self):
        cases = {case for case, _ in self.table}
        self.assertEqual(cases, {'SSSS', 'F_y0', 'F_x0', 'FF_y', 'FF_x',
                                 'FF_adj', 'FFF_x', 'FFF_y', 'FFFF'})
        self.assertEqual(len(self.table), 81)
        self.assertEqual(sum(len(x) for x in self.table.values()), 1701)
        for (case, quantity), rows in self.table.items():
            self.assertEqual(len(rows), 21, (case, quantity))
            self.assertEqual([round(r['ratio_ly_lx'], 2) for r in rows],
                             [round(1 + .05*i, 2) for i in range(21)])
            self.assertTrue(all(r['table_value'] > 0 and
                                math.isfinite(r['table_value']) for r in rows))

    def test_square_simple_moment(self):
        factor, meta = lookup('SSSS', 1, 'mx_mid', table=self.table)
        self.assertFalse(meta['interpolated'])
        self.assertAlmostEqual(factor, 1 / 27.2)
        self.assertAlmostEqual(calculate('SSSS', 4, 4, 10, 33e6, .18,
                                         table=self.table)['responses']['mx_mid']['value'],
                               160 / 27.2)

    def test_square_fixed_negative(self):
        factor, meta = lookup('FFFF', 1, 'mx_fixed_mid', table=self.table)
        self.assertEqual(meta['sign'], '-')
        self.assertAlmostEqual(factor, -1 / 19.4)

    def test_corner_force_consistency_square(self):
        f, _ = lookup('SSSS', 1, 'corner_force', table=self.table)
        t, _ = lookup('SSSS', 1, 'mxy_corner', table=self.table)
        self.assertAlmostEqual(f, 2 * t)

    def test_new_005_nodes_are_exact_not_interpolated(self):
        anchors = [
            ('F_y0', 1.05, 'mx_mid', 36.5),
            ('F_x0', 1.10, 'mx_mid', 27.3),
            ('FF_y', 1.65, 'my_fixed_mid', 9.05),
            ('FF_x', 1.55, 'mx_fixed_mid', 12.1),
            ('FF_adj', 1.85, 'mx_mid', 18.7),
            ('FFF_x', 1.95, 'my_mid', 95.4),
            ('FFF_y', 1.90, 'my_max', 91.7),
        ]
        for case, ratio, quantity, denominator in anchors:
            with self.subTest(case=case, ratio=ratio, quantity=quantity):
                _, meta = lookup(case, ratio, quantity, table=self.table)
                self.assertFalse(meta['interpolated'])
                self.assertAlmostEqual(meta['table_value'], denominator)

    def test_interpolate_only_between_005_knots(self):
        mid, meta = lookup('F_y0', 1.125, 'mx_mid', table=self.table)
        self.assertTrue(meta['interpolated'])
        self.assertAlmostEqual(mid, .5/31.9 + .5/28.3)
        by_denominator, _ = lookup('F_y0', 1.125, 'mx_mid',
                                   table=self.table, interpolate='denominator')
        self.assertAlmostEqual(by_denominator, 1/((31.9+28.3)/2))

    def test_interpolation_not_extrapolation(self):
        with self.assertRaises(ValueError):
            lookup('SSSS', 2.01, 'mx_mid', table=self.table)
        with self.assertRaises(ValueError):
            lookup('FF_adj', 0.99, 'mx_mid', table=self.table)


if __name__ == '__main__':
    unittest.main()
