"""Tests for the pure/geometry-utility functions in common/__init__.py.

None of these touch VASP -- they exercise math helpers (get_line,
get_parab), force/geometry helpers (fmax, get_closest, pos_swap,
reindex_atoms, get_nearest_neighbors, match_cell), and the INCAR text
editor (param_set).
"""
import numpy as np
import pytest

import common


def test_get_line_recovers_slope_intercept():
    x = [0.0, 1.0, 2.0, 3.0]
    m, b = 2.0, 1.0
    y = [m * xi + b for xi in x]

    xax, yax, a, bb = common.get_line(x, y)
    assert a == pytest.approx(m, abs=1e-8)
    assert bb == pytest.approx(b, abs=1e-8)
    # xax should span beyond [min(x), max(x)] by `extra`
    assert xax.min() < min(x)
    assert xax.max() > max(x)


def test_get_line_return_mae_is_zero_for_exact_line():
    x = [0.0, 1.0, 2.0, 3.0, 4.0]
    y = [3.0 * xi - 2.0 for xi in x]
    xax, yax, a, b, mae = common.get_line(x, y, return_mae=True)
    assert a == pytest.approx(3.0, abs=1e-8)
    assert b == pytest.approx(-2.0, abs=1e-8)
    assert mae == pytest.approx(0.0, abs=1e-8)


def test_get_parab_endpoints_match():
    x0, x1 = 0.0, 1.0
    y0, y1 = 0.0, 1.5

    # side 0: left half of the parabola (zero slope at x1)
    xax, yax = common.get_parab(x0, x1, y0, y1, side=0)
    assert xax[0] == pytest.approx(x0)
    assert xax[-1] == pytest.approx(x1)
    assert yax[0] == pytest.approx(y0, abs=1e-6)
    assert yax[-1] == pytest.approx(y1, abs=1e-6)

    # side 1: right half of the parabola (zero slope at x0)
    xax, yax = common.get_parab(x0, x1, y0, y1, side=1)
    assert yax[0] == pytest.approx(y0, abs=1e-6)
    assert yax[-1] == pytest.approx(y1, abs=1e-6)


def _make_atoms_with_forces(forces, fixed_indices=None):
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.constraints import FixAtoms

    n = len(forces)
    positions = [(i * 2.0, 0.0, 0.0) for i in range(n)]
    atoms = Atoms('H%d' % n, positions=positions, cell=[10, 10, 10], pbc=False)
    if fixed_indices:
        atoms.set_constraint(FixAtoms(indices=fixed_indices))
    calc = SinglePointCalculator(atoms, energy=0.0, forces=np.array(forces))
    atoms.calc = calc
    return atoms


def test_fmax_ignores_fixed_atoms():
    # atom 0 has the largest force but is fixed; atom 1's force should win
    forces = [
        [10.0, 0.0, 0.0],  # fixed -- must be excluded
        [0.0, 3.0, 4.0],   # magnitude 5.0
        [0.0, 0.0, 1.0],   # magnitude 1.0
    ]
    atoms = _make_atoms_with_forces(forces, fixed_indices=[0])
    assert common.fmax(atoms) == pytest.approx(5.0, abs=1e-8)


def test_fmax_no_constraints():
    forces = [[3.0, 4.0, 0.0], [1.0, 0.0, 0.0]]
    atoms = _make_atoms_with_forces(forces)
    assert common.fmax(atoms) == pytest.approx(5.0, abs=1e-8)


def test_get_closest():
    from ase import Atoms

    ref = Atoms('H2', positions=[(0, 0, 0), (5, 0, 0)], cell=[20, 20, 20], pbc=False)
    atoms = Atoms('H2', positions=[(0.1, 0, 0), (4.9, 0, 0)], cell=[20, 20, 20], pbc=False)

    # atoms[0] (near 0,0,0) should be closest to ref[0]
    assert common.get_closest(ref, atoms, 0, mic=False) == 0
    # atoms closest to ref[1] (near 5,0,0) should be atoms[1]
    assert common.get_closest(ref, atoms, 1, mic=False) == 1


def test_pos_swap():
    from ase import Atoms

    atoms = Atoms('H2', positions=[(1, 2, 3), (4, 5, 6)], cell=[20, 20, 20], pbc=False)
    common.pos_swap(atoms, 0, 1)
    assert atoms[0].position == pytest.approx([4, 5, 6])
    assert atoms[1].position == pytest.approx([1, 2, 3])


def test_reindex_atoms_recovers_shuffled_order():
    from ase import Atoms

    ref = Atoms('H3', positions=[(0, 0, 0), (2, 0, 0), (4, 0, 0)], cell=[20, 20, 20], pbc=False)
    # shuffled copy: index 0 <-> index 2 swapped relative to ref
    shuffled = Atoms('H3', positions=[(4, 0, 0), (2, 0, 0), (0, 0, 0)], cell=[20, 20, 20], pbc=False)

    reindexed = common.reindex_atoms(ref, shuffled)
    for i in range(len(ref)):
        assert reindexed[i].position == pytest.approx(ref[i].position, abs=1e-6)


def test_get_nearest_neighbors_dimer_within_and_outside_cutoff():
    from ase import Atoms
    from ase.data import covalent_radii, atomic_numbers

    r_h = covalent_radii[atomic_numbers['H']]
    scale = 1.1

    # within cutoff: bond length just under scale*(r_h+r_h)
    d_in = scale * (2 * r_h) * 0.9
    atoms_close = Atoms('H2', positions=[(0, 0, 0), (d_in, 0, 0)], cell=[20, 20, 20], pbc=False)
    nn = common.get_nearest_neighbors(atoms_close, scale=scale)
    assert 1 in nn[0]
    assert 0 in nn[1]

    # outside cutoff: bond length well beyond scale*(r_h+r_h)
    d_out = scale * (2 * r_h) * 3.0
    atoms_far = Atoms('H2', positions=[(0, 0, 0), (d_out, 0, 0)], cell=[20, 20, 20], pbc=False)
    nn_far = common.get_nearest_neighbors(atoms_far, scale=scale)
    assert nn_far[0] == []
    assert nn_far[1] == []


def test_match_cell_anchor_atom():
    from ase import Atoms

    ref = Atoms('H2', positions=[(0, 0, 5), (0, 0, 6)], cell=[10, 10, 20], pbc=False)
    change = Atoms('H2', positions=[(0, 0, 8), (0, 0, 9)], cell=[5, 5, 5], pbc=False)

    out = common.match_cell(ref, change, lower_vac=0.0, anchor_atom=0)
    # cell should now match ref
    assert np.allclose(np.array(out.cell), np.array(ref.cell))
    # anchor atom (index 0) should now share z with ref's anchor atom
    assert out[0].z == pytest.approx(ref[0].z, abs=1e-8)
    # the rigid z-shift between the two atoms should be preserved
    assert (out[1].z - out[0].z) == pytest.approx(change[1].z - change[0].z, abs=1e-8)


def test_match_cell_lower_vac():
    from ase import Atoms

    ref = Atoms('H1', positions=[(0, 0, 0)], cell=[10, 10, 20], pbc=False)
    change = Atoms('H2', positions=[(0, 0, 3), (0, 0, 8)], cell=[5, 5, 5], pbc=False)

    out = common.match_cell(ref, change, lower_vac=1.0)
    zs = sorted(atom.z for atom in out)
    assert zs[0] == pytest.approx(1.0, abs=1e-8)


def test_param_set_modifies_existing_tag(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'INCAR').write_text(
        "ENCUT = 400\nISMEAR = 0\n"
    )
    common.param_set('ENCUT', 520)
    lines = (tmp_path / 'INCAR').read_text().splitlines()
    encut_lines = [l for l in lines if 'ENCUT' in l]
    assert len(encut_lines) == 1
    assert '520' in encut_lines[0]
    # untouched tag survives
    assert any('ISMEAR' in l for l in lines)


def test_param_set_appends_missing_tag(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'INCAR').write_text("ENCUT = 400\n")
    common.param_set('NSW', 300)
    lines = (tmp_path / 'INCAR').read_text().splitlines()
    nsw_lines = [l for l in lines if 'NSW' in l]
    assert len(nsw_lines) == 1
    assert '300' in nsw_lines[0]


def test_param_set_del_removes_tag(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'INCAR').write_text("ENCUT = 400\nNSW = 300\n")
    common.param_set('NSW', 'del')
    lines = (tmp_path / 'INCAR').read_text().splitlines()
    assert not any('NSW' in l for l in lines)
    assert any('ENCUT' in l for l in lines)
