"""Template: constant-potential geometry relaxation.

Copy this into a job directory and edit the potential, kpts, and any VASP
tags below for your system. Calls common.const_U_relax so the NELECT/
Newton-step logic and checkpointing live in one place instead of being
re-implemented per job directory.
"""

from common import const_U_relax
from ase.calculators.vasp import Vasp
from ase.io import read
import os

kpts = (4, 4, 1)
pot_des = -1.413  # desired potential, V vs SHE
ediffg = 0.05

# resume from a previous CONTCAR if present, otherwise start from POSCAR
try:
    atoms = read('CONTCAR')
except Exception:
    atoms = read('POSCAR')

calc = Vasp(txt='-',
            encut=500,
            xc='PBE',
            gga='RP',
            kpts=kpts,
            ncore=16,
            kpar=1,
            gamma=True,  # Gamma-centered (defaults to Monkhorst-Pack)
            ismear=0,
            algo='Normal',
            nelm=100,
            sigma=0.05,
            ibrion=1,
            ediffg=-ediffg,  # forces
            ediff=1e-4,
            prec='Accurate',
            nsw=0,
            lcharg=False,
            lsol=True,
            lambda_d_k=3.0,
            tau=0.0)

# const_U_relax handles set_pot + geometry optimization + checkpointing.
# Default optimizer='vasp' uses VASP's internal IBRION=1 optimizer, same as
# above. To instead hand geometry steps to the VTST optimizer compiled into
# VASP (IBRION=3/POTIM=0/IOPT), pass optimizer='vtst' (iopt=2 for CG,
# iopt=7 for FIRE), e.g.:
#   const_U_relax(atoms, calc, pot_des, ediffg=ediffg, optimizer='vtst', iopt=7)
const_U_relax(atoms, calc, pot_des, ediffg=ediffg)
