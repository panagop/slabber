"""Manual transcription of *uniformly distributed full load* entries from
J. Hahn's Tafel 17--22 (printed pp. 373--379). One series per image/table row.
NOT claimed to be Czerny's complete original 1996 tables.
Read the source and THREE_EDGE_SOURCES.md before using these data.
"""
from __future__ import annotations
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = 'https://fgg-web.fgg.uni-lj.si/KMLK/Sebastjan/masivni%20objekti/Hahnove%20tabele.pdf'
R14 = '1.50 1.40 1.30 1.20 1.10 1.00 0.90 0.80 0.70 0.60 0.50 0.40 0.30 0.25'
R11 = '2.00 1.80 1.60 1.40 1.20 1.00 0.80 0.60 0.50 0.40 0.25'
# tuple: printed symbol, quantity type, sign, value kind, scale, values in ratio order, note
# "denominator" => response = scale / table_value, multiplier => scale * table_value.
# `K = p*lx*ly`. All moment coefficients are Hahn's K / k format.
CASES = {
'H17_SSS': (17,374,'S,FREE,S,S',R14, {
 'mx_free_mid': ('m_xr','moment','+','denominator','K', '12.6 11.9 11.3 10.7 10.2 9.8 9.4 9.1 9.1 9.2 9.8 11.0 13.7 16.2',''),
 'mx_center': ('m_xm','moment','+','denominator','K', '15.3 14.9 14.5 14.1 13.8 13.7 13.6 13.8 14.2 15.2 17.0 20.2 26.3 31.5',''),
 'my_center': ('m_ym','moment','+','denominator','K', '62.4 58.4 54.2 50.0 45.9 41.7 37.1 33.2 29.9 27.4 25.9 26.3 29.7 33.7',''),
 'mxy_rear_corner': ('m_xy2','moment','±','denominator','K', '22.3 20.6 19.3 17.9 16.7 15.4 14.1 12.9 11.8 10.8 10.1 9.4 8.8 8.6','Rear two-supported-edge corner; R2=2*m_xy2'),
 'mxy_free_corner': ('m_xy1','moment','±','denominator','K', '412 300 220 161 118 86.5 63.6 47.0 35.0 26.3 20.2 15.8 12.8 11.6','Free/supported corner; R1=2*m_xy1'),
 'deflection_ref': ('omega_r','deflection','+','denominator','K_lx2_over_Ed3','9.10 8.70 8.35 8.05 7.80 7.60 7.45 7.35 7.35 7.40 7.65 8.25 9.90 11.60','Hahn: w_r = K*lx^2/(omega_bar*E*h^3); D = omega_bar*E*h^3'),
 'reaction_x_total': ('nu_x','reaction','+','multiplier','K','0.45 0.45 0.44 0.43 0.42 0.41 0.39 0.37 0.34 0.31 0.28 0.22 0.16 0.13','K_x=nu_x*K, resultant for each x-side edge, see source sketch'),
 'reaction_y_total': ('nu_y','reaction','+','multiplier','K','0.28 0.30 0.32 0.34 0.36 0.40 0.44 0.49 0.54 0.59 0.64 0.72 0.80 0.84','K_y=nu_y*K, resultant at rear edge, see source sketch'),
}),
'H18_REAR_FIXED': (18,375,'S,FREE,S,F',R14,{
 'mx_free_mid': ('m_xr','moment','+','denominator','K','13.1 12.5 12.1 11.7 11.5 11.4 11.5 12.0 13.0 15.2 19.4 29.4 60.2 105',''),
 'mx_center': ('m_xm','moment','+','denominator','K','18.1 18.1 18.1 18.3 18.8 19.7 21.0 23.3 27.0 34.2 48.0 79.0 174 293',''),
 'my_positive_adjacent': ("m_y'm",'moment','+','denominator','K','84 77 70 64 59 55 52 54 57 63 72 85 107 124','Positive moment near fixed rear edge'),
 'my_rear_fixed': ('m_ey','moment','-','denominator','K','12.1 11.3 10.5 9.8 9.1 8.5 7.9 7.4 7.1 6.8 6.8 7.1 8.1 9.0',''),
 'mxy_free_corner': ('m_xy1','moment','±','denominator','K','262 195 146 110 84 64 48 40 33 29 26 26 30 35','R=2*m_xy1, magnitude'),
 'negative_moment_reach': ('beta','fraction','+','multiplier','1','0.20 0.22 0.25 0.28 0.31 0.34 0.37 0.40 0.43 0.46 0.49 0.53 0.57 0.60','Approximate distance beta*ly of the negative my zone from the rear support'),
}),
'H19_ONE_SIDE_FIXED': (19,376,'F,FREE,S,S',R11,{
 'mx_free_mid': ('m_xr','moment','+','denominator','K','27.8 26.3 23.2 20.4 17.7 15.3 12.9 11.7 11.4 11.6 15.4',''),
 'mx_center': ('m_xm','moment','+','denominator','K','30.7 27.5 24.8 22.5 20.5 18.8 17.5 17.7 18.5 20.6 31.3',''),
 'my_center_or_max': ('m_y','moment','+','denominator','K','117 101 86 71 58 46 35 28 29 30 35','For eps > approximately 0.8, this represents my_max near rear edge, not my at center'),
 'mx_fixed_free_corner': ('m_er','moment','-','denominator','K','14.6 13.3 12.0 10.6 9.2 7.9 6.7 5.8 5.0 4.5 4.1',''),
 'mx_fixed_side_mid': ('m_em','moment','-','denominator','K','17.0 15.4 13.9 12.5 11.2 10.1 9.2 8.8 9.1 9.8 11.7',''),
}),
'H20_BOTH_SIDES_FIXED': (20,377,'F,FREE,F,S',R11,{
 'mx_free_mid': ('m_xr','moment','+','denominator','K','48.1 43.3 38.4 33.4 28.3 23.5 19.5 15.9 14.8 13.8 16.6',''),
 'mx_center': ('m_xm','moment','+','denominator','K','49.0 44.6 40.2 36.4 32.5 29.0 25.4 22.7 22.1 23.3 29.0',''),
 'my_center_or_max': ('m_y','moment','+','denominator','K','170 149 127 108 90 72 52 44 39 36 37','For eps > approximately 0.6, this represents my_max near rear edge, not my at center'),
 'mx_fixed_free_corner': ('m_er','moment','-','denominator','K','24.1 21.6 19.2 16.9 14.3 11.9 9.5 7.1 6.1 5.2 4.5',''),
 'mx_fixed_side_mid': ('m_em','moment','-','denominator','K','24.2 21.9 19.7 17.4 15.7 13.8 12.1 10.7 10.1 9.7 9.8',''),
}),
'H21_ADJ_FIXED': (21,378,'F,FREE,S,F',R11,{
 'mx_free_mid': ('m_xr','moment','+','denominator','K','30.1 26.7 23.7 21.0 18.6 16.4 15.1 16.0 19.1 25.2 54.0',''),
 'mx_center': ('m_xm','moment','+','denominator','K','33.2 30.3 27.3 25.7 24.8 24.9 26.8 33.3 41.4 55.5 112',''),
 'my_center': ('m_ym','moment','+','denominator','K','201 195 137 125 106 83 67 71 80 103 690','Check rapid printed variation near eps 0.25 and 1.8 against physical copy before design'),
 'mx_fixed_free_corner': ('m_er','moment','-','denominator','K','14.7 13.3 11.9 10.5 9.1 7.8 6.7 6.0 5.7 5.8 6.7',''),
 'mx_fixed_side_mid': ('m_em','moment','-','denominator','K','17.9 16.3 14.8 13.5 12.6 12.1 12.3 13.6 15.1 17.4 26.1',''),
 'my_rear_fixed': ('m_ey','moment','-','denominator','K','25.3 22.9 20.5 18.0 15.4 12.8 10.3 8.4 8.0 7.9 9.1',''),
}),
'H22_ALL_FIXED': (22,379,'F,FREE,F,F',R14,{
 'mx_free_mid': ('m_xr','moment','+','denominator','K','35.8 33.4 31.0 28.6 26.4 24.3 22.4 20.9 19.9 19.8 21.3 26.8 46.4 77.0',''),
 'mx_center': ('m_xm','moment','+','denominator','K','39.8 38.3 37.0 35.8 34.9 34.3 34.0 34.3 35.6 38.6 45.6 63.6 126 228',''),
 'my_center': ('m_ym','moment','+','denominator','K','163 152 141 130 119 109 99.5 91.0 83.4 80.0 83.4 108 208 417',''),
 'mx_fixed_free_corner': ('m_er','moment','-','denominator','K','17.8 16.6 15.3 14.1 12.8 11.6 10.4 9.3 8.2 7.4 6.8 6.8 7.6 8.6',''),
 'mx_fixed_side_mid': ('m_em','moment','-','denominator','K','18.7 17.8 17.0 16.2 15.6 15.0 14.5 14.3 14.2 14.7 15.8 18.1 23.0 27.2',''),
 'my_rear_fixed': ('m_ey','moment','-','denominator','K','26.4 24.6 22.8 21.1 19.3 17.6 15.8 14.2 12.6 11.1 9.8 9.0 9.0 9.6',''),
}),
}

COLUMNS = ['case_id','ratio_ly_lx','quantity','original_symbol','measure_kind','sign','table_value','value_kind','scale_basis','source_table','source_page','source_pdf_page','source_url','notes']
def build(out=ROOT / 'czerny_three_edge_coefficients.csv'):
    rows=[]
    for case_id, (no, page, edges, ratios_text, data) in CASES.items():
        ratios=[float(s) for s in ratios_text.split()]
        assert len(ratios)==len(set(ratios)) and all(b<a for a,b in zip(ratios,ratios[1:]))
        for quantity, (symbol,kind,sign,valkind,scale,values,notes) in data.items():
            values=[float(s) for s in values.split()]
            assert len(ratios)==len(values), (case_id,quantity,len(ratios),len(values))
            assert all(v>0 for v in values), (case_id,quantity)
            for ratio,value in zip(ratios,values):
                rows.append(dict(case_id=case_id,ratio_ly_lx=f'{ratio:.2f}',quantity=quantity,original_symbol=symbol,measure_kind=kind,sign=sign,table_value=f'{value:.9g}',value_kind=valkind,scale_basis=scale,source_table=f'Hahn Tafel {no}',source_page=page,source_pdf_page=page-361,source_url=SOURCE,notes=notes))
    rows.sort(key=lambda r: (r['case_id'],r['quantity'],float(r['ratio_ly_lx'])))
    with open(out,'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=COLUMNS)
        w.writeheader();w.writerows(rows)
    print(f'Created {out}: {len(rows)} coefficient entries across {len(CASES)} cases')
    return rows

if __name__=='__main__':
    build()
