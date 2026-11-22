"""Smoke/regression checks against printed Hahn 1970 Tafeln 17-22 (uniform load only).
These tests do NOT independently verify every manually transcribed cell.
"""
import math
import unittest

from czerny import calculate, load_table, load_all_tables, lookup, THREE_EDGE_CSV


class HahnThreeEdgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hahn = load_table(THREE_EDGE_CSV)
        cls.both = load_all_tables()

    def test_all_six_cases_and_456_original_nodes(self):
        cases={case for case,_ in self.hahn}
        self.assertEqual(cases,{
            'H17_SSS','H18_REAR_FIXED','H19_ONE_SIDE_FIXED',
            'H20_BOTH_SIDES_FIXED','H21_ADJ_FIXED','H22_ALL_FIXED',
        })
        self.assertEqual(sum(len(v) for v in self.hahn.values()),456)
        self.assertEqual(len(self.hahn),36)
        for (case,name),rows in self.hahn.items():
            expected=14 if case in ('H17_SSS','H18_REAR_FIXED','H22_ALL_FIXED') else 11
            self.assertEqual(len(rows),expected,(case,name))
            self.assertTrue(all(r['table_value']>0 and math.isfinite(r['table_value']) for r in rows))
            self.assertTrue(all(r['ratio_ly_lx']<s['ratio_ly_lx'] for r,s in zip(rows,rows[1:])))

    def test_four_edge_kept_separate_and_all_loadable(self):
        self.assertEqual(len(self.both),81+36)
        self.assertEqual(sum(len(r) for r in self.both.values()),1701+456)
        self.assertEqual(len(load_table()),81)

    def test_hahn_table_17_printed_example_ratio_0_6(self):
        # Source printed p. 374: lx=2.5, ly=1.5, p=9.5 -> K=35.625kN.
        K=9.5*2.5*1.5
        r=calculate('H17_SSS',2.5,1.5,9.5,33e6,.18,table=self.both)
        self.assertAlmostEqual(r['total_uniform_load_K'],K)
        self.assertAlmostEqual(r['responses']['mx_free_mid']['value'],K/9.2)
        self.assertAlmostEqual(r['responses']['mx_center']['value'],K/15.2)
        self.assertAlmostEqual(r['responses']['my_center']['value'],K/27.4)
        self.assertAlmostEqual(r['responses']['mxy_rear_corner']['value'],K/10.8)
        self.assertAlmostEqual(r['responses']['mxy_free_corner']['value'],K/26.3)
        self.assertTrue(r['responses']['mxy_free_corner']['value_is_magnitude_only'])
        self.assertAlmostEqual(r['responses']['reaction_x_total']['value'],K*.31)
        self.assertAlmostEqual(r['responses']['reaction_y_total']['value'],K*.59)
        self.assertAlmostEqual(r['responses']['deflection_ref']['value'],
                               K*2.5**2/(7.4*33e6*.18**3))

    def test_hahn_table_18_beta(self):
        r=calculate('H18_REAR_FIXED',3.,1.8,10.,33e6,.18,table=self.both)
        # epsilon .6; beta = 0.46 => length 0.828m
        self.assertAlmostEqual(r['responses']['negative_moment_reach']['value'],.46)
        self.assertAlmostEqual(r['responses']['negative_moment_reach']['extent_from_rear_m'],.828)
        self.assertAlmostEqual(r['responses']['my_rear_fixed']['value'],-54/6.8)

    def test_fixed_cases_printed_anchors(self):
        anchors=[
            ('H19_ONE_SIDE_FIXED',.8,'mx_fixed_side_mid',9.2),
            ('H20_BOTH_SIDES_FIXED',1.2,'my_center_or_max',90),
            ('H21_ADJ_FIXED',1.0,'my_rear_fixed',12.8),
            ('H22_ALL_FIXED',.6,'my_rear_fixed',11.1),
        ]
        for case,ratio,qty,k in anchors:
            with self.subTest(case=case,ratio=ratio,qty=qty):
                factor,meta=lookup(case,ratio,qty,table=self.both)
                self.assertFalse(meta['interpolated'])
                self.assertEqual(meta['scale_basis'],'K')
                self.assertEqual(meta['table_value'],k)
                self.assertAlmostEqual(factor,-1/k if qty.startswith('mx_fixed') or qty=='my_rear_fixed' else 1/k)

    def test_non_tabulated_ratio_interpolation_and_bounds(self):
        factor,meta=lookup('H22_ALL_FIXED',.65,'mx_free_mid',table=self.both)
        self.assertTrue(meta['interpolated'])
        self.assertAlmostEqual(factor,.5/19.8+.5/19.9)
        self.assertEqual(meta['source_table'],'Hahn Tafel 22')
        with self.assertRaises(ValueError):
            lookup('H22_ALL_FIXED',2.0,'mx_free_mid',table=self.both)
        with self.assertRaises(ValueError):
            lookup('H20_BOTH_SIDES_FIXED',.2,'mx_free_mid',table=self.both)

    def test_results_only_for_values_the_source_provides(self):
        r=calculate('H22_ALL_FIXED',2.5,1.5,7.5,33e6,.18,table=self.both)
        self.assertAlmostEqual(r['responses']['mx_free_mid']['value'],7.5*2.5*1.5/19.8)
        self.assertNotIn('qx_mid',r['responses'])
        self.assertNotIn('deflection_ref',r['responses'])
        self.assertEqual(len(r['responses']),6)
        self.assertEqual(r['responses']['mx_fixed_side_mid']['unit'],'kNm/m')

    def test_zero_and_negative_load_signs(self):
        x=calculate('H19_ONE_SIDE_FIXED',2.5,1.5,0,33e6,.18,table=self.both)
        self.assertTrue(all(v['value']==0 for v in x['responses'].values()))
        y=calculate('H19_ONE_SIDE_FIXED',2.5,1.5,-10,33e6,.18,table=self.both)
        self.assertGreater(y['responses']['mx_fixed_free_corner']['value'],0)

if __name__=='__main__':
    unittest.main()
