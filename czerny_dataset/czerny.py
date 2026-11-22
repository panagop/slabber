"""Published slab-coefficient lookup: Czerny four-edge tables and Hahn three-edge tables.

* Four edge: Czerny, Beton-Kalender 1996 pp. 281-289; CSV `czerny_coefficients.csv`.
  moment/corner p*lx**2 / k; shear p*lx / k; deflection c*p*lx**4/(E*h**3).
* Three edge: Hahn, Tafeln 17-22 pp. 373-379; CSV
  `czerny_three_edge_coefficients.csv`; moments K/k with K=p*lx*ly.
  Hahn's uniform-load table does NOT give a full shear/deflection set for all cases.

Research/transcription aid only, NOT a validated design tool. The three-edge
values are from HAHN (secondary compilation), not a recovered full Czerny 1996
original. Both table series presume a linear elastic isotropic constant-thickness
plate with the source's Poisson-ratio assumptions. No extrapolation allowed.

Units: lx, ly, thickness [m]; load [kN/m²]; E [kN/m²]. Moments [kNm/m],
reactions and corner forces [kN], shear [kN/m], displacement [m].

Coordinates for 3-edge cases: x = 0..lx horizontally, y = 0..ly vertically;
y=0 is the FREE edge, y=ly is rear supported edge. The lx dimension is not
necessarily the shorter span for three-edge cases, unlike four-edge cases.
Three-edge edge ordering is (x=0, y=0, x=lx, y=ly). `F` means FIXED;
`FREE` means free/unrestrained.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DEFAULT_CSV = ROOT / 'czerny_coefficients.csv'
THREE_EDGE_CSV = ROOT / 'czerny_three_edge_coefficients.csv'


def load_table(path=DEFAULT_CSV):
    """Load one long-form CSV, keyed by (case_id, quantity)."""
    result = {}
    with open(path, newline='', encoding='utf-8-sig') as file:
        for row in csv.DictReader(file):
            row['ratio_ly_lx'] = float(row['ratio_ly_lx'])
            row['table_value'] = float(row['table_value'])
            if row['table_value'] <= 0 or not math.isfinite(row['table_value']):
                raise ValueError(f'Invalid coefficient: {row}')
            key = (row['case_id'], row['quantity'])
            result.setdefault(key, []).append(row)
    for key, rows in result.items():
        rows.sort(key=lambda x: x['ratio_ly_lx'])
        if any(a['ratio_ly_lx'] >= b['ratio_ly_lx']
               for a, b in zip(rows, rows[1:])):
            raise ValueError(f'Duplicate aspect ratio for {key}')
    return result


def load_all_tables(*paths):
    """Load both 4-edge Czerny and 3-edge Hahn (or selected CSV paths)."""
    if not paths:
        paths = (DEFAULT_CSV, THREE_EDGE_CSV)
    result = {}
    for path in paths:
        for key, values in load_table(path).items():
            if key in result:
                raise ValueError(f'Duplicate key from CSVs: {key}')
            result[key] = values
    return result


def _basis(row):
    """Return the unit scaling convention of each original table."""
    if row.get('scale_basis'):
        return row['scale_basis']
    kind = row['measure_kind']  # legacy Czerny schema
    if kind in ('moment', 'corner_force'):
        return 'p_lx2'
    if kind == 'shear':
        return 'p_lx'
    if kind == 'deflection':
        return 'p_lx4_over_Ed3'
    raise ValueError(f'Unknown source scaling for {row}')


def lookup(case_id, ratio, quantity, *, table=None, interpolate='coefficient'):
    """Return (signed physical multiplier, metadata) for a tabulated result.

    Printed denominator k => signed factor 1/k. Printed direct multiplier c
    => signed factor c. For signs marked ± the returned number is the positive
    MAGNITUDE, without positional/directional sign information.

    Default interpolation is linear in physical response factor (1/k for
    denominator entries). Use `interpolate='denominator'` for linear in k.
    Interpolation is NOT source-certified and does not extrapolate.
    """
    if not math.isfinite(ratio):
        raise ValueError('Aspect ratio must be finite')
    if interpolate not in ('coefficient', 'denominator'):
        raise ValueError('interpolate must be coefficient or denominator')
    data = load_all_tables() if table is None else table
    samples = data.get((case_id, quantity))
    if samples is None:
        raise KeyError(f'Unknown case/quantity: {case_id} / {quantity}')
    xs = [d['ratio_ly_lx'] for d in samples]
    if ratio < xs[0] - 1e-10 or ratio > xs[-1] + 1e-10:
        raise ValueError(f'{case_id}/{quantity}: ratio outside [{xs[0]}, {xs[-1]}]')
    exactly = next((r for r in samples if abs(r['ratio_ly_lx']-ratio)<1e-10), None)
    if exactly is not None:
        row=exactly
        value=row['table_value']
        interpolated=False
    else:
        lo,hi=next((a,b) for a,b in zip(samples,samples[1:])
                   if a['ratio_ly_lx']<ratio<b['ratio_ly_lx'])
        if (lo['value_kind'],lo['measure_kind'],lo['sign'],_basis(lo)) != (hi['value_kind'],hi['measure_kind'],hi['sign'],_basis(hi)):
            raise ValueError('Inconsistent coefficient conventions in table')
        t=(ratio-lo['ratio_ly_lx'])/(hi['ratio_ly_lx']-lo['ratio_ly_lx'])
        if lo['value_kind']=='denominator' and interpolate=='coefficient':
            factor=(1-t)/lo['table_value'] + t/hi['table_value']
            value=1/factor
        else:
            value=(1-t)*lo['table_value'] + t*hi['table_value']
        row=lo
        interpolated=True
    if row['value_kind'] == 'denominator':
        factor=1/value
    elif row['value_kind']=='multiplier':
        factor=value
    else:
        raise ValueError(f'Unknown value_kind: {row["value_kind"]}')
    sign=row['sign']
    return (-factor if sign=='-' else factor), {
        'case_id':case_id,
        'quantity':quantity,
        'ratio_ly_lx':ratio,
        'original_symbol':row['original_symbol'],
        'sign':sign,
        'measure_kind':row['measure_kind'],
        'scale_basis':_basis(row),
        'table_value':value,
        'table_value_kind':row['value_kind'],
        'interpolated':interpolated,
        'source_page':int(row['source_page']),
        'source_table':row.get('source_table','Czerny 1996 four-edge'),
        'source_url':row.get('source_url',''),
        'notes':row.get('notes',''),
    }


def calculate(case_id, lx, ly, load, E, thickness, *, table=None,
              interpolate='coefficient'):
    """Calculate ALL *available*, not all theoretically possible, table results.

    E must be >0 even if selected case has no deflection entry to preserve
    the uniform API. Output keeps directional limitations and source provenance.
    """
    if not all(math.isfinite(float(v)) for v in (lx,ly,load,E,thickness)):
        raise ValueError('All inputs must be finite')
    if min(lx,ly,E,thickness)<=0:
        raise ValueError('Spans, E and thickness must be >0')
    ratio=ly/lx
    data=load_all_tables() if table is None else table
    names=sorted(quantity for case,quantity in data if case==case_id)
    if not names:
        raise KeyError(f'Unsupported case: {case_id}')
    K=load*lx*ly
    result={
        'case_id':case_id,
        'ratio_ly_lx':ratio,
        'total_uniform_load_K':K,
        'three_edges_supported':case_id.startswith('H'),
        'nu_assumed':0,
        'warning':'Elastic uncracked plate coefficients only; do not use for final structural design. ± signs denote magnitudes only. See documentation.',
        'responses':{},
    }
    scales={
        'p_lx2': load*lx**2,
        'p_lx':load*lx,
        'p_lx4_over_Ed3':load*lx**4/(E*thickness**3),
        'K':K,
        'K_lx2_over_Ed3':K*lx**2/(E*thickness**3),
        '1':1,
    }
    for name in names:
        multiplier, meta=lookup(case_id,ratio,name,table=data,interpolate=interpolate)
        kind=meta['measure_kind']
        units={'moment':'kNm/m','corner_force':'kN','reaction':'kN',
               'shear':'kN/m','deflection':'m','fraction':'dimensionless'}
        if kind not in units:
            raise ValueError(f'Unknown measure_kind: {kind}')
        basis=meta['scale_basis']
        if basis not in scales:
            raise ValueError(f'Unknown scale_basis: {basis}')
        response={
            'value':scales[basis]*multiplier,
            'unit':units[kind],
            'sign':meta['sign'],
            'value_is_magnitude_only':meta['sign']=='±',
            'interpolated':meta['interpolated'],
            'source_table':meta['source_table'],
            'source_page':meta['source_page'],
            'printed_table_value':meta['table_value'],
            'printed_value_kind':meta['table_value_kind'],
            'scale_basis':basis,
        }
        if meta['notes']:
            response['notes']=meta['notes']
        if name=='negative_moment_reach':
            response['extent_from_rear_m']=response['value']*ly
        result['responses'][name]=response
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',default='SSSS',help='SSSS, FFFF, H17_SSS, ... H22_ALL_FIXED')
    ap.add_argument('--lx',type=float,default=4.0)
    ap.add_argument('--ly',type=float,default=6.0)
    ap.add_argument('--load',type=float,default=10.0,help='uniform p [kN/m^2]')
    ap.add_argument('--E',type=float,default=33e6,help='elastic modulus [kN/m^2]')
    ap.add_argument('--thickness',type=float,default=.18)
    ap.add_argument('--interpolate',choices=['coefficient','denominator'],default='coefficient')
    ap.add_argument('--csv',type=Path,nargs='*',help='optional one or multiple CSV files')
    args=ap.parse_args()
    table=load_all_tables(*args.csv) if args.csv else load_all_tables()
    print(json.dumps(calculate(args.case,args.lx,args.ly,args.load,args.E,args.thickness,
                               table=table,interpolate=args.interpolate),
                     indent=2,ensure_ascii=False))

if __name__=='__main__':
    main()
