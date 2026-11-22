"""Compute all 21 SSSS Czerny ratios at one mesh density with scipy FD."""
import csv
from pathlib import Path
from finite_difference_reference import solve_ssss
from czerny import lookup

ROOT=Path(__file__).resolve().parent
rows=[]
for i in range(21):
    ratio=round(1+i*.05,2)
    nx=40;ny=40+2*i
    v=solve_ssss(nx=nx,ny=ny,b=4*ratio)
    kx=lookup('SSSS',ratio,'mx_mid')[1]['table_value']
    ky=lookup('SSSS',ratio,'my_max')[1]['table_value']
    c=lookup('SSSS',ratio,'deflection_mid')[1]['table_value']
    v.update(kx_czerny=kx,ky_czerny=ky,c_czerny=c,
             kx_error_percent=100*(v['kx']/kx-1),
             ky_error_percent=100*(v['ky']/ky-1),
             c_error_percent=100*(v['deflection_multiplier']/c-1))
    rows.append(v)
with open(ROOT/'reference_ssss_fd_all_ratios.csv','w',newline='') as out:
    w=csv.DictWriter(out,list(rows[0]));w.writeheader();w.writerows(rows)
for field in ['kx_error_percent','ky_error_percent','c_error_percent']:
    mx=max(rows,key=lambda v:abs(v[field]))
    print(field,'max_abs=',round(mx[field],3),'at ratio',mx['ratio_ly_lx'])
print('Rows:',len(rows))
