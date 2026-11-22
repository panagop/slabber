"""PyNiteFEA 3.2.x reproduction experiment for Czerny/Hahn rectangular plates.

Install locally: python -m pip install "PyNiteFEA==3.2.1"
Run from this directory:
  python pynite_czerny_benchmark.py --case SSSS --ratio 1.5 --meshes 7 11 15

REPRODUCIBILITY STATUS: the code has been syntax-checked but this environment
could not install PyNiteFEA, so this module has NOT been solver-executed here.
All numbers in reference_ssss_fd.csv come from a separate SciPy finite-difference
benchmark, NOT from PyNite.

Plate lies in XY plane; downward load is negative local/global Z.
Edge ordering: left/x=0, bottom/y=0, right/x=a, rear/y=b.
S = simply supported (w=0; rotations free), F = fixed/clamped (w=0,
  both rotations blocked), FREE = traction-free (out-of-plane DOFs free).
For the *flat, unloaded-in-plane* pure-bending benchmark only, DX, DY and RZ
are restrained at EVERY node to eliminate irrelevant membrane/drilling modes.
This is NOT a model of in-plane slab structural action.

Output moment extrema are samples at element CENTERS; this is not identical
to the exact continuous maximum, and near-edge/corner extrema are sensitive
to mesh and sampling. Source moment and shear sign conventions may differ.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def edge_codes(case_id: str):
    with (ROOT/'coverage.csv').open(newline='', encoding='utf-8') as stream:
        for row in csv.DictReader(stream):
            if row['case_id']==case_id:
                codes=row['edges_x0_y0_xa_yb'].split(',')
                if len(codes)!=4 or not set(codes).issubset({'S','F','FREE'}):
                    raise ValueError(f'Unsupported edge convention: {codes}')
                return codes
    raise KeyError(f'No such case in coverage.csv: {case_id}')


def nearest_odd(n):
    n=max(3,round(n))
    return n if n%2 else n+1


def solve_pynite(case_id='SSSS',a=4.0,b=6.0,p=10.0,E=33e6,h=0.18,
                 nu=0.0,nx=11):
    """Return a PyNite result report from a flat rectangular Rect-element mesh."""
    try:
        from Pynite import FEModel3D
    except ImportError as exc:
        raise RuntimeError('PyNiteFEA is required: pip install PyNiteFEA==3.2.1') from exc
    if nx<3 or nx%2!=1:
        raise ValueError('Use an odd nx >=3 so the slab center is inside an element')
    if not all(math.isfinite(v) for v in (a,b,p,E,h,nu)):
        raise ValueError('All numerical inputs must be finite')
    if min(a,b,p,E,h)<=0 or not (-1<nu<0.5):
        raise ValueError('Invalid plate geometry/material/load')

    left,bottom,right,rear=edge_codes(case_id)
    ny=nearest_odd(nx*b/a)
    dx,dy=a/nx,b/ny
    m=FEModel3D()
    m.add_material('Concrete',E,G=E/(2*(1+nu)),nu=nu,rho=0.0)
    node_id=lambda i,j:f'N_{i}_{j}'

    # Generate a regular, consistently wound rectangular element mesh.
    for j in range(ny+1):
        for i in range(nx+1):
            m.add_node(node_id(i,j),i*dx,j*dy,0.0)

    for j in range(ny):
        for i in range(nx):
            name=f'P_{i}_{j}'
            m.add_plate(name, node_id(i,j),node_id(i+1,j),
                        node_id(i+1,j+1),node_id(i,j+1),h,'Concrete')
            m.add_plate_surface_pressure(name, -p, case='P')

    for j in range(ny+1):
        for i in range(nx+1):
            flags=[]
            if i==0: flags.append(left)
            if j==0: flags.append(bottom)
            if i==nx: flags.append(right)
            if j==ny: flags.append(rear)
            support_z=any(s in ('S','F') for s in flags)
            clamp=any(s=='F' for s in flags)
            m.def_support(node_id(i,j),support_DX=True,support_DY=True,
                          support_DZ=support_z,support_RX=clamp,
                          support_RY=clamp,support_RZ=True)

    m.add_load_combo('P',{'P':1.0})
    m.analyze_linear(log=False,check_statics=True)

    # Midspan is the CENTER of one element because nx,ny are odd.
    ci,cj=nx//2,ny//2
    mid_plate=m.plates[f'P_{ci}_{cj}']
    mx_mid,my_mid,mxy_mid=(float(x) for x in
        mid_plate.moment(x=dx/2,y=dy/2,local=True,combo_name='P').flatten())
    qx_mid,qy_mid=(float(x) for x in
        mid_plate.shear(x=dx/2,y=dy/2,local=True,combo_name='P').flatten())

    # 4 nodal displacements averaged at the geometric center; this interpolates
    # the center deflection and converges on mesh refinement.
    wmid=sum(m.nodes[node_id(i,j)].DZ['P']
             for i in (ci,ci+1) for j in (cj,cj+1))/4

    sampled=[]
    for j in range(ny):
        for i in range(nx):
            plate=m.plates[f'P_{i}_{j}']
            mx,my,mxy=(float(x) for x in
                plate.moment(x=dx/2,y=dy/2,local=True,combo_name='P').flatten())
            sampled.append(dict(x=(i+0.5)*dx,y=(j+0.5)*dy,
                                Mx=mx,My=my,Mxy=mxy))

    # Sampling ALL element centers catches off-center bending extrema; signs
    # are not interpreted until compared with the source's sign convention.
    ext={}
    for name in ('Mx','My','Mxy'):
        most_pos=max(sampled,key=lambda row:row[name])
        most_neg=min(sampled,key=lambda row:row[name])
        most_abs=max(sampled,key=lambda row:abs(row[name]))
        ext[name]={
            'max':{k:most_pos[k] for k in ('x','y',name)},
            'min':{k:most_neg[k] for k in ('x','y',name)},
            'max_abs':{k:most_abs[k] for k in ('x','y',name)},
        }

    rz=sum(float(n.RxnFZ['P']) for n in m.nodes.values())
    sssp=case_id=='SSSS'
    table={}
    if sssp and abs(nu)<1e-12:
        # Czerny denominator normalization m=p*lx^2/k.
        from czerny import lookup
        base=p*a*a
        for quantity,value in [('mx_mid',mx_mid),
                               ('my_max',ext['My']['max_abs']['My'])]:
            c,meta=lookup(case_id,b/a,quantity)
            table[quantity]={
                'source_k':meta['table_value'],
                'fem_k':base/abs(value),
                'moment':value,
                'relative_error_percent':100*(abs(value)/(base*c)-1),
            }
        _,meta=lookup(case_id,b/a,'deflection_mid')
        c_fem=abs(wmid)*E*h**3/(p*a**4)
        table['deflection_mid']={
            'source_multiplier':meta['table_value'],
            'fem_multiplier':c_fem,
            'relative_error_percent':100*(c_fem/meta['table_value']-1),
        }

    return {
        'solver':'PyNiteFEA 3.2.x (Rect plate elements)',
        'case_id':case_id,'edges_x0_y0_xa_yb':[left,bottom,right,rear],
        'lx':a,'ly':b,'ratio_ly_lx':b/a,'uniform_load':p,
        'E':E,'h':h,'nu':nu,'nx':nx,'ny':ny,'elements':nx*ny,
        'midspan':{'Mx':mx_mid,'My':my_mid,'Mxy':mxy_mid,
                   'Qx':qx_mid,'Qy':qy_mid,'w':wmid},
        'sampled_extrema':ext,
        'reaction_sum_kN':rz,
        'applied_total_load_kN':p*a*b,
        'absolute_reaction_balance_percent':100*(abs(rz)/(p*a*b)-1),
        'czerny_comparison':table,
        'notes':'Do not use extrema near point support/corners without a separate convergence check.',
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',default='SSSS')
    ap.add_argument('--lx',type=float,default=4.)
    ap.add_argument('--ratio',type=float,default=1.5)
    ap.add_argument('--load',type=float,default=10.)
    ap.add_argument('--E',type=float,default=33e6)
    ap.add_argument('--h',type=float,default=0.18)
    ap.add_argument('--nu',type=float,default=0)
    ap.add_argument('--meshes',nargs='+',type=int,default=[7,11,15])
    ap.add_argument('--output',type=Path,default=ROOT/'pynite_results.json')
    args=ap.parse_args()
    rows=[]
    for nx in args.meshes:
        res=solve_pynite(case_id=args.case,a=args.lx,b=args.lx*args.ratio,
                         p=args.load,E=args.E,h=args.h,nu=args.nu,nx=nx)
        rows.append(res)
        print('case=',args.case,'nx=',res['nx'],'ny=',res['ny'],
              'reaction balance %=',res['absolute_reaction_balance_percent'])
        print('  midspan:',res['midspan'])
        print('  Czerny comparison:',res['czerny_comparison'])
    args.output.write_text(json.dumps(rows,indent=2),encoding='utf-8')
    print('Saved:',args.output)

if __name__=='__main__':
    main()
