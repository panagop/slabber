"""Independent finite-difference benchmark for all-four-edge simply-supported slabs.

Only the SSSS case with nu=0 is covered. This is NOT PyNite and not a
replacement for the complete Czerny tables. It solves D*(Laplacian^2 w)=p
with Navier boundary conditions w=0 and Laplacian(w)=0 on straight edges.

Units: a,b,h m; p kN/m2; E kN/m2.  Outputs m, kNm/m, etc.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import factorized


def solve_ssss(a=4., b=4., p=10., E=33e6, thickness=0.18,
               nx=32, ny=None, nu=0.0):
    if nu != 0:
        raise ValueError('This Navier FD implementation currently validates nu=0 only')
    if nx < 4 or nx % 2:
        raise ValueError('nx must be even and >=4, to sample exactly at midspan')
    if ny is None:
        ny=round(nx*b/a)
    if ny < 4 or ny % 2:
        raise ValueError('ny must be even and >=4, to sample exactly at midspan')
    hx,hy=a/nx,b/ny
    nxi,nyi=nx-1,ny-1
    Dx=diags([np.ones(nxi-1),-2*np.ones(nxi),np.ones(nxi-1)],[-1,0,1],format='csc')/hx**2
    Dy=diags([np.ones(nyi-1),-2*np.ones(nyi),np.ones(nyi-1)],[-1,0,1],format='csc')/hy**2
    L=kron(eye(nyi,format='csc'),Dx,format='csc') + kron(Dy,eye(nxi,format='csc'),format='csc')
    D=E*thickness**3/(12*(1-nu**2))
    linear=factorized(L)
    U=linear(np.full(nxi*nyi, p/D))
    W=linear(U).reshape((nyi,nxi))
    ix=nx//2-1;iy=ny//2-1
    # Center derivatives and global extrema from the interior nodal field.
    padded=np.pad(W,((1,1),(1,1)),constant_values=0.0)
    Mx=-D*(padded[1:-1,2:]-2*W+padded[1:-1,:-2])/hx**2
    My=-D*(padded[2:,1:-1]-2*W+padded[:-2,1:-1])/hy**2
    mx=Mx[iy,ix]
    my=My[iy,ix]
    jmax,imax=np.unravel_index(My.argmax(),My.shape)
    mymax=My[jmax,imax]
    w=W[iy,ix]
    return dict(nx=nx,ny=ny,element_size_x=hx,element_size_y=hy,
                lx=a,ly=b,ratio_ly_lx=b/a,load=p,nu=nu,
                Mx_mid=mx, My_mid=my, My_max=mymax,
                My_max_x=(imax+1)*hx,My_max_y=(jmax+1)*hy,
                w_mid=w,
                kx=p*a*a/mx,ky=p*a*a/mymax,
                deflection_multiplier=w*E*thickness**3/(p*a**4))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ratios',nargs='+',type=float,default=[1.0,1.5,2.0])
    parser.add_argument('--refinements',nargs='+',type=int,default=[8,16,32,64])
    parser.add_argument('--lx',type=float,default=4)
    parser.add_argument('--load',type=float,default=10)
    parser.add_argument('--E',type=float,default=33e6)
    parser.add_argument('--thickness',type=float,default=0.18)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('reference_ssss_fd.csv'))
    args=parser.parse_args()
    from czerny import lookup
    rows=[]
    for ratio in args.ratios:
        base_kx=lookup('SSSS',ratio,'mx_mid')[1]['table_value']
        base_ky=lookup('SSSS',ratio,'my_max')[1]['table_value']
        base_w=lookup('SSSS',ratio,'deflection_mid')[1]['table_value']
        for nx in args.refinements:
            v=solve_ssss(a=args.lx,b=args.lx*ratio,p=args.load,E=args.E,thickness=args.thickness,nx=nx)
            v.update(kx_czerny=base_kx,ky_czerny=base_ky,
                     cx_czerny=base_w,
                     mx_diff_percent=100*(v['Mx_mid']/(args.load*args.lx**2/base_kx)-1),
                     my_diff_percent=100*(v['My_max']/(args.load*args.lx**2/base_ky)-1),
                     w_diff_percent=100*(v['deflection_multiplier']/base_w-1))
            rows.append(v)
    with args.output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(f'Saved {len(rows)} comparisons in {args.output}')
    for row in rows:
        print(f"a/b={row['ratio_ly_lx']:.2f}, nx={row['nx']:3d}, "
              f"Mx={row['Mx_mid']:.5f}, My={row['My_mid']:.5f}, "
              f"kx={row['kx']:.4f}, ky={row['ky']:.4f}, "
              f"mx err={row['mx_diff_percent']:+.2f}%, my err={row['my_diff_percent']:+.2f}%, "
              f"defl err={row['w_diff_percent']:+.2f}%")

if __name__=='__main__':
    main()
