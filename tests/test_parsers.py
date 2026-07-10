"""Fixture-based tests for the file-parsing functions in common/__init__.py.

All fixture files below are fabricated-minimal -- just enough content to
exercise the parser's actual split()/substring logic (read the parser
first: it grabs `split()[2]` off any line containing 'fermi', a
' = '-delimited value off any line containing 'FERMI_SHIFT', etc). They
are not real VASP output and can be swapped for real recorded files later.

No test here invokes VASP; everything runs against files written directly
in tmp_path.
"""
import pytest

import common


def test_get_wf_implicit_vaspsolpp_isol_branch(tmp_path):
    # fabricated-minimal OUTCAR: one 'fermi' line (parser wants
    # split()[2] to be the fermi value) plus an ISOL line, which routes
    # get_wf_implicit down the VASPsol++ branch (no FERMI_SHIFT lookup).
    d = tmp_path / 'isol_job'
    d.mkdir()
    (d / 'OUTCAR').write_text(
        " E-fermi :  -2.1235     XC(G=0):  -1.2345  alpha+bet : -0.5678\n"
        " ISOL  =    1\n"
    )
    wf = common.get_wf_implicit(str(d))
    assert wf == pytest.approx(2.1235)


def test_get_wf_implicit_vaspsol_fermi_shift_branch(tmp_path):
    # fabricated-minimal OUTCAR (no ISOL line) plus an opt.log carrying
    # FERMI_SHIFT, which routes get_wf_implicit down the classic VASPsol
    # branch: -1*(fermi+shift).
    d = tmp_path / 'sol_job'
    d.mkdir()
    (d / 'OUTCAR').write_text(
        " E-fermi :  -3.5000     XC(G=0):  -1.0000  alpha+bet : -0.2000\n"
    )
    (d / 'opt.log').write_text("FERMI_SHIFT = 0.0500\n")
    wf = common.get_wf_implicit(str(d))
    # -1*(fermi + shift) = -1*(-3.5 + 0.05) = 3.45
    assert wf == pytest.approx(3.45)


def test_get_wf_explicit(tmp_path):
    # this also regression-tests the Item A fix: get_wf_explicit uses
    # os.path.isfile but previously never imported os.
    d = tmp_path / 'wf_job'
    d.mkdir()
    (d / 'wf.out').write_text("1.2345\n")
    wf = common.get_wf_explicit(str(d))
    assert wf == pytest.approx(1.2345)


def test_get_omega_no_outcar_no_vasprun_path(tmp_path):
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointCalculator

    d = tmp_path / 'omega_job'
    d.mkdir()

    atoms = Atoms('H2', positions=[(0, 0, 0), (0, 0, 0.74)], cell=[10, 10, 10], pbc=False)

    # POSCAR so get_n0 has a geometry file to read (n0 = 2 * ZVAL(H) = 2.0)
    atoms.write(str(d / 'POSCAR'))

    # lastimage.traj carries the energy get_omega falls back to reading
    # once there's no OUTCAR/vasprun.xml
    calc = SinglePointCalculator(atoms, energy=-10.0)
    atoms.calc = calc
    atoms.write(str(d / 'lastimage.traj'))

    (d / 'nel.txt').write_text("10.0\n")
    (d / 'fermi.txt').write_text("-4.0\n")

    omega = common.get_omega(str(d))
    # e - q*fermi, q = nel - n0 = 10.0 - 2.0 = 8.0
    # omega = -10.0 - 8.0*(-4.0) = 22.0
    assert omega == pytest.approx(22.0)
