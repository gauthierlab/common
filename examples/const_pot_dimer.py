"""Template: constant-potential Dimer transition-state search.

Copy this into a job directory and edit the potential, kpts, and any VASP
tags below for your system. Calls common.const_U_dimer so the NELECT/
Newton-step logic and checkpointing live in one place instead of being
re-implemented per job directory.

Ideally you should already have a MODECAR file in the job directory --
convergence is poor without a good initial guess.
"""

from common import const_U_dimer
from ase.calculators.vasp import Vasp
from ase.io import read

kpts = (4, 4, 1)
pot_des = -1.413  # desired potential, V vs SHE
ediffg = 0.02

# CENTCAR is the correct restart geometry for the dimer method;
# fall back to POSCAR for a fresh start
try:
    atoms = read('CENTCAR')
except Exception:
    atoms = read('POSCAR')

calc = Vasp(txt='-',
            encut=500,
            xc='PBE',
            gga='RP',
            kpts=kpts,
            ncore=16,
            kpar=3,
            gamma=True,  # Gamma-centered (defaults to Monkhorst-Pack)
            ismear=0,
            algo='Normal',
            nelm=100,
            sigma=0.05,
            ediffg=-ediffg,  # forces
            ediff=1e-8,
            prec='Accurate',
            nsw=0,
            lcharg=False,
            lsol=True,
            lambda_d_k=3.0,
            tau=0.0,
            lasph=True,
            ibrion=3,
            potim=0,
            ichain=2,
            iopt=2)

# const_U_dimer handles set_pot + Dimer optimization + checkpointing
# (CENTCAR/MODECAR housekeeping included). iopt=2 is CG (default); 7 is
# FIRE, which can be more robust than the L-BFGS VTST default (iopt=1).
const_U_dimer(atoms, calc, pot_des, ediff=1e-8, ediffg=ediffg, iopt=2)
