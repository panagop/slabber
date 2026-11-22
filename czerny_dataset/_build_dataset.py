import csv, json, pathlib, zipfile
root=pathlib.Path(__file__).parent
SOURCE='https://www.klancek.si/sites/default/files/datoteke/files/tabele.pdf'
ratios21=[round(1+i*.05,2) for i in range(21)]
ratios5=[1.,1.25,1.5,1.75,2.]
CASES={
'SSSS': (281, 'S,S,S,S', ratios21, {
'mx_mid':('m_xzm','moment','+', '27.2 24.5 22.4 20.7 19.1 17.8 16.8 15.8 15.0 14.3 13.7 13.2 12.7 12.3 11.9 11.5 11.3 11.0 10.8 10.6 10.4'),
'my_max':('m_ymax','moment','+', '27.2 27.5 27.9 28.4 29.1 29.9 30.9 31.8 32.8 33.8 34.7 35.4 36.1 36.7 37.3 37.9 38.5 38.9 39.4 39.8 40.3'),
'mxy_corner':('m_xye','moment','±', '21.6 20.6 19.7 19.0 18.4 17.9 17.5 17.1 16.8 16.5 16.3 16.1 15.9 15.7 15.6 15.5 15.4 15.3 15.3 15.2 15.1'),
'corner_force':('R_e','corner_force','±', '10.8 10.3 9.85 9.5 9.2 8.95 8.75 8.55 8.4 8.25 8.15 8.05 7.95 7.85 7.8 7.75 7.7 7.65 7.65 7.6 7.55'),
'qx_mid':('q_xrm','shear','±', '2.96 2.87 2.78 2.71 2.64 2.58 2.52 2.47 2.43 2.39 2.36 2.33 2.30 2.27 2.25 2.23 2.21 2.19 2.18 2.16 2.15'),
'qx_support_mid':('qbar_xrm','shear','+', '2.19 2.15 2.11 2.07 2.04 2.02 2.00 1.98 1.97 1.96 1.95 1.94 1.93 1.92 1.92 1.92 1.92 1.92 1.92 1.92 1.92'),
'qy_mid':('q_yrm','shear','±', '2.96 2.92 2.89 2.86 2.84 2.82 2.80 2.78 2.76 2.75 2.75 2.74 2.73 2.73 2.73 2.72 2.72 2.71 2.71 2.70 2.70'),
'qy_support_mid':('qbar_yrm','shear','+', '2.19 2.14 2.09 2.05 2.02 1.99 1.96 1.94 1.92 1.90 1.89 1.88 1.87 1.86 1.85 1.84 1.83 1.82 1.82 1.82 1.82'),
'deflection_mid':('f_m','deflection','+', '.0487 .0536 .0584 .0631 .0678 .0728 .0767 .0809 .0850 .0890 .0927 .0963 .0997 .1029 .1060 .1093 .1118 .1145 .1169 .1195 .1215'),
}),
'F_y0': (282,'S,F,S,S',ratios5,{
'mx_mid':('m_xzm','moment','+','41.2 23.4 16.6 13.2 11.4'),
'my_fixed_mid':('m_yerm','moment','-','11.9 9.8 8.9 8.45 8.2'),
'my_max':('m_ymax','moment','+','29.4 29.2 32.3 35.6 38.8'),
'mxy_corner':('m_xye','moment','±','26.2 20.1 17.4 16.1 15.4'),
'corner_force':('R_e','corner_force','±','13.1 10.0 8.7 8.0 7.7'),
'qx_mid':('q_xrm','shear','±','3.54 2.92 2.57 2.36 2.23'),
'qx_support_mid':('qbar_xrm','shear','+','2.59 2.23 2.07 1.98 1.95'),
'qy_fixed_mid':('q_yerm','shear','+','1.72 1.53 1.44 1.39 1.37'),
'qy_mid':('q_yrm','shear','-','3.18 2.92 2.80 2.74 2.71'),
'qy_support_mid':('qbar_yrm','shear','+','2.47 2.12 1.96 1.88 1.84'),
'deflection_mid':('f_m','deflection','+','.0334 .0557 .0773 .0962 .1112'),
}),
'F_x0': (283,'F,S,S,S',ratios5,{
'mx_fixed_mid':('m_xerm','moment','-','11.9 9.9 9.0 8.5 8.3'),
'mx_mid':('m_xzm','moment','+','31.4 23.4 19.8 18.0 17.1'),
'my_max':('m_ymax','moment','+','41.2 50.3 55.6 58.2 59.2'),
'mxy_corner':('m_xye','moment','±','26.2 23.7 22.8 22.4 22.4'),
'corner_force':('R_e','corner_force','±','13.1 11.8 11.4 11.2 11.2'),
'qx_fixed_mid':('q_xerm','shear','+','1.72 1.62 1.58 1.57 1.57'),
'qx_mid':('q_xrm','shear','-','3.18 2.90 2.77 2.70 2.67'),
'qx_support_mid':('qbar_xrm','shear','+','2.47 2.40 2.42 2.45 2.50'),
'qy_mid':('q_yrm','shear','±','3.53 3.53 3.54 3.55 3.56'),
'qy_support_mid':('qbar_yrm','shear','+','2.59 2.50 2.48 2.47 2.47'),
'deflection_mid':('f_m','deflection','+','.0334 .0438 .0510 .0556 .0585'),
}),
'FF_y':(284,'S,F,S,F',ratios5,{
'mx_mid':('m_xzm','moment','+','63.3 31.5 20.4 15.2 12.5'),
'my_fixed_mid':('m_yerm','moment','-','14.3 11.1 9.5 8.8 8.4'),
'my_max':('m_ymax','moment','+','35.1 31.3 32.1 35.9 42.4'),
'qx_mid':('q_xrm','shear','±','4.09 3.27 2.78 2.49 2.32'),
'qx_support_mid':('qbar_xrm','shear','+','2.95 2.44 2.18 2.05 1.98'),
'qy_fixed_mid':('q_yerm','shear','±','1.94 1.65 1.50 1.42 1.38'),
'deflection_mid':('f_m','deflection','+','.0230 .0425 .0639 .0840 .1013'),
}),
'FF_x':(285,'F,S,F,S',ratios5,{
'mx_fixed_mid':('m_xerm','moment','-','14.3 12.7 12.2 12.0 12.0'),
'mx_mid':('m_xzm','moment','+','35.1 28.5 25.8 24.5 24.1'),
'my_max':('m_ymax','moment','+','61.7 72.8 75.8 77.0 77.0'),
'qx_fixed_mid':('q_xerm','shear','±','1.94 1.90 1.91 1.93 1.95'),
'qy_mid':('q_yrm','shear','±','4.09 4.14 4.17 4.19 4.19'),
'qy_support_mid':('qbar_yrm','shear','+','2.95 2.95 2.96 2.96 2.97'),
'deflection_mid':('f_m','deflection','+','.0230 .0275 .0297 .0308 .0313'),
}),
'FF_adj':(286,'F,F,S,S',ratios5,{
'mx_fixed_mid':('m_xerm','moment','-','14.3 11.1 9.6 8.8 8.4'),
'mx_mid':('m_xzm','moment','+','42.7 28.0 22.2 19.5 17.9'),
'my_fixed_mid':('m_yerm','moment','-','14.3 12.9 12.4 12.2 12.2'),
'my_max':('m_ymax','moment','+','40.2 45.6 53.0 57.0 60.2'),
'qx_fixed_mid':('q_xerm','shear','+','1.96 1.72 1.64 1.61 1.58'),
'qx_mid':('q_xrm','shear','-','3.64 3.14 2.89 2.76 2.71'),
'qx_support_mid':('qbar_xrm','shear','+','2.76 2.53 2.46 2.47 2.49'),
'qy_fixed_mid':('q_yerm','shear','+','1.96 1.86 1.83 1.82 1.82'),
'qy_mid':('q_yrm','shear','-','3.64 3.54 3.54 3.54 3.55'),
'qy_support_mid':('qbar_yrm','shear','+','2.76 2.56 2.49 2.48 2.47'),
'deflection_mid':('f_m','deflection','+','.0252 .0369 .0459 .0519 .0562'),
}),
'FFF_x':(287,'F,F,S,F',ratios5,{
'mx_fixed_mid':('m_xerm','moment','-','18.3 12.7 10.6 9.5 8.8'),
'mx_mid':('m_xzm','moment','+','59.5 34.2 25.2 21.1 18.7'),
'my_fixed_min':('m_yermin','moment','-','16.2 13.5 12.7 12.4 12.3'),
'my_mid':('m_ym','moment','+','44.1 45.8 55.0 75.0 101.0'),
'qx_fixed_mid':('q_xerm','shear','+','2.21 1.85 1.71 1.66 1.61'),
'qx_mid':('q_xrm','shear','-','4.11 3.40 3.03 2.85 2.75'),
'qx_support_mid':('qbar_xrm','shear','+','3.04 2.65 2.51 2.48 2.47'),
'qy_fixed_mid':('q_yerm','shear','±','2.11 1.91 1.85 1.83 1.82'),
'deflection_mid':('f_m','deflection','+','.0188 .0308 .0409 .0484 .0539'),
}),
'FFF_y':(288,'F,F,F,S',ratios5,{
'mx_fixed_min':('m_xermin','moment','-','16.2 13.5 12.5 12.1 12.0'),
'mx_mid':('m_xzm','moment','+','44.1 32.3 27.6 25.3 24.5'),
'my_fixed_mid':('m_yerm','moment','-','18.3 17.5 17.5 17.5 17.5'),
'my_max':('m_ymax','moment','+','55.9 67.7 75.2 84.6 97.0'),
'qx_fixed_mid':('q_xerm','shear','±','2.11 1.94 1.94 1.94 1.95'),
'qy_fixed_mid':('q_yerm','shear','+','2.21 2.10 2.19 2.20 2.21'),
'qy_mid':('q_yrm','shear','-','4.11 4.12 4.15 4.15 4.15'),
'qy_support_mid':('qbar_yrm','shear','+','3.04 2.79 2.78 2.78 2.78'),
'deflection_mid':('f_m','deflection','+','.0188 .0245 .0280 .0301 .0309'),
}),
'FFFF':(289,'F,F,F,F',ratios21,{
'mx_fixed_mid':('m_xerm','moment','-','19.4 18.2 17.1 16.3 15.5 14.9 14.5 14.0 13.7 13.4 13.2 13.0 12.8 12.7 12.5 12.4 12.3 12.2 12.1 12.0 12.0'),
'mx_mid':('m_xzm','moment','+','56.8 50.6 46.1 42.4 39.4 37.0 34.8 33.3 31.9 30.6 29.6 28.8 28.1 27.5 26.9 26.4 26.0 25.7 25.4 25.2 25.0'),
'my_fixed_mid':('m_yerm','moment','-','19.4 18.8 18.4 18.1 17.9 17.7 17.6 17.5 17.5 17.5 17.5 17.5 17.5 17.5 17.5 17.5 17.5 17.5 17.5 17.5 17.5'),
'my_max':('m_ymax','moment','+','56.8 58.2 60.3 62.6 65.8 69.4 73.6 78.4 83.4 89.4 93.5 96.1 98.1 99.9 101.3 102.4 103.3 104.0 104.6 104.9 105.0'),
'qx_fixed_mid':('q_xerm','shear','±','2.24 2.17 2.10 2.05 2.01 1.98 1.96 1.94 1.92 1.92 1.92 1.92 1.92 1.92 1.91 1.91 1.91 1.91 1.91 1.91 1.91'),
'qy_fixed_mid':('q_yerm','shear','±','2.24 2.20 2.16 2.14 2.12 2.11 2.10 2.09 2.09 2.10 2.12 2.12 2.12 2.12 2.12 2.12 2.12 2.12 2.13 2.13 2.13'),
'deflection_mid':('f_m','deflection','+','.0152 .0167 .0181 .0195 .0207 .0219 .0230 .0240 .0248 .0257 .0264 .0271 .0277 .0282 .0288 .0291 .0294 .0297 .0300 .0302 .0304'),
}),
}
# Complete all seven previously sampled cases from the printed 0.05 grid.
# Data are transcribed in _full_grid.py, separately from original five-point
# seed samples, so both sources can be cross-checked in code review.
from _full_grid import FULL
for _case, _fields in FULL.items():
    _page, _edges, _old_ratios, _original = CASES[_case]
    assert set(_fields) == set(_original), (_case, 'missing field')
    _expanded = {}
    for _key, _raw in _fields.items():
        _values = [float(_v) for _v in _raw.split()]
        assert len(_values) == 21, (_case, _key, len(_values))
        _symbol, _kind, _sign, _old = _original[_key]
        _expanded[_key] = (_symbol, _kind, _sign, _raw)
    CASES[_case] = (_page, _edges, ratios21, _expanded)

cols=['case_id','ratio_ly_lx','quantity','original_symbol','measure_kind','sign','table_value','value_kind','source_page','source_url']
rows=[]; coverage=[]
for case,(page,edge,ratios,values) in CASES.items():
    for field,(symbol,kind,sign,vals) in values.items():
        num=list(map(float,vals.split()))
        assert len(num)==len(ratios),(case,field,len(num),len(ratios))
        for ratio,v in zip(ratios,num):
            rows.append(dict(case_id=case,ratio_ly_lx=f'{ratio:.2f}',quantity=field,original_symbol=symbol,measure_kind=kind,sign=sign,table_value=f'{v:.8g}',value_kind='multiplier' if kind=='deflection' else 'denominator',source_page=page,source_url=SOURCE))
    coverage.append(dict(case_id=case,edges_x0_y0_xa_yb=edge,source_page=page,aspect_ratios=len(ratios),ratio_min=min(ratios),ratio_max=max(ratios),n_parameters=len(values),completeness='full 0.05 grid' if len(ratios)==21 else 'sampled 0.25 grid',source_url=SOURCE))
with (root/'czerny_coefficients.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(rows)
with (root/'coverage.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=coverage[0].keys()); w.writeheader();w.writerows(coverage)
print('Rows:',len(rows),'Cases:',len(CASES),'ratios:',{k:len(v[2]) for k,v in CASES.items()})
