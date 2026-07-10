"""Template: constant-potential fixed-bond-length (FBL) optimization.

Copy this into a job directory and edit the potential, atom indices, and
kpts below for your system. Calls common.const_U_FBL so the NELECT/Newton-
step logic and checkpointing live in one place instead of being
re-implemented per job directory.

FBL is an alternative to the Dimer method when your reaction pathway is
roughly one dimensional -- it fixes the distance between two atoms
(ind1, ind2) via ase.constraints.FixBondLength while everything else
relaxes.
"""

from common import const_U_FBL
from ase.calculators.vasp import Vasp
from ase.io import read

kpts = (4, 4, 1)
pot_des = -1.413  # desired potential, V vs SHE
ediffg = 0.05

o_ind = 36
h_ind = 37
z_cutoff = 12.2  # fix atoms below this z, in addition to the FBL pair

# restart order: CONTCAR (mid-relaxation), then init.traj (freshly set up
# by e.g. common.get_interp), then POSCAR
try:
    atoms = read('CONTCAR')
except Exception:
    try:
        atoms = read('init.traj')
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

# const_U_FBL handles set_pot + FixBondLength/FixAtoms constraints +
# ASE-side optimization + checkpointing. Default ase_optimizer='bfgs'
# matches previous behavior; if BFGS struggles to converge, try
# ase_optimizer='fire' (same algorithm as VTST IOPT=7, but enforced on
# the ASE side so the bond constraint stays intact), e.g.:
#   const_U_FBL(atoms, calc, pot_des, o_ind, h_ind, z_cutoff=z_cutoff,
#               ediffg=ediffg, ase_optimizer='fire')
const_U_FBL(atoms, calc, pot_des, o_ind, h_ind, z_cutoff=z_cutoff, ediffg=ediffg)
