"""Convergence of the supercell deflection efficiency (LPCVD family, p = 0.6 um, h = 1.6 um, 1075 nm).

Usage: python scripts/supercell_convergence.py > results/supercell_convergence/log.txt
"""
import time, numpy as np
from metatrap import materials
from metatrap.rcwa.library import Library
from metatrap.rcwa import lookup as L
from metatrap.rcwa.supercell import deflection_efficiency
lib=Library.load('data/library_main.h5'); lam=1.075
n_in=float(materials.FUSED_SILICA.n(lam)); n_p=complex(materials.SIN_LPCVD.nk(lam))
p,h=0.6,1.6
w,t,_=lib.curve('sin_lpcvd',p,h,lam); lk=L.build(w,t,p,h,n_in,1.0,L.FabLimits())
for nc in [4,8]:
    for tpc,ng in [(80,48),(140,48),(200,64),(140,96)]:
        t0=time.time(); r=deflection_efficiency(lk,nc,p,h,lam,n_p**2,n_in**2,n_grid=ng,terms_per_cell=tpc,polarisation='x')
        print(f'N={nc} terms/cell={tpc} grid={ng}: eta={r["eta_rcwa"]:.4f} R={r["R"]:.4f} eta0={r["eta_zero"]:.4f} ({time.time()-t0:.0f}s)', flush=True)
# uniform-array check: N cells all the same width -> should reproduce single-cell T0, R
from metatrap.rcwa.solver import solve_grating_orders, solve_pillars
from metatrap.rcwa.supercell import supercell_eps
from metatrap.rcwa.cell import PillarCell
wmid=0.33
eps=supercell_eps([wmid]*4,p,n_p**2,1.0,48)
o,e,R=solve_grating_orders(eps,4*p,p,h,lam,n_in**2,1.0,n_terms=80*4)
s=solve_pillars(PillarCell(p,h),[wmid],lam,n_p**2,n_in**2,n_terms=197,n_grid=256)
print('uniform 4-cell: T0 %.4f R %.4f | single cell: T0 %.4f R %.4f'%(dict(zip(map(tuple,o),e))[(0,0)],R,s.T0[0],s.R[0]))
