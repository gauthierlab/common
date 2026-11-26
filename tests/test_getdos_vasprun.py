"""get_dos: DOSCAR path vs vasprun.xml fallback on a synthetic 2-atom LORBIT 11 file pair."""
import os
import pickle

import numpy as np
import pytest

import common

ORBS = ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2']
NEDOS, EF = 6, -5.1234


def _data(ispin, seed=0):
    rng = np.random.default_rng(seed)
    e = np.round(np.linspace(-10.0, 2.0, NEDOS), 3)
    tot = np.round(rng.random((NEDOS, ispin)) * 3, 4)
    pd = np.round(rng.random((2, ispin, 9, NEDOS)), 4)
    return e, tot, np.cumsum(tot, axis=0), pd


def _write_pair(d, ispin, total_only=False):
    e, tot, tin, pd = _data(ispin)
    from ase import Atoms
    from ase.io import write
    write(os.path.join(d, 'POSCAR'), Atoms('AuO', positions=[[0, 0, 0], [0, 0, 2]], cell=[5, 5, 10]),
          format='vasp')
    L = ['   2   2   1   0', ' 0.1E+03 0.1E-08 0.1E-08 0.1E-08 0.5E-15', ' 1.0E-004', '  CAR ', ' x',
         '  %12.8f %12.8f %d %12.8f 1.0' % (e[-1], e[0], NEDOS, EF)]
    for i in range(NEDOS):
        L.append(' %10.3f' % e[i] + ''.join(' %11.4E' % v for v in tot[i]) + ''.join(' %11.4E' % v for v in tin[i]))
    if not total_only:
        for a in range(2):
            L.append(L[5])
            for i in range(NEDOS):
                cols = [pd[a, s, o, i] for o in range(9) for s in range(ispin)]
                L.append(' %10.3f' % e[i] + ''.join(' %11.4E' % v for v in cols))
    open(os.path.join(d, 'DOSCAR'), 'w').write('\n'.join(L) + '\n')

    X = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<modeling>',
         '<atominfo><atoms> 2 </atoms><array name="atoms"><set>',
         '<rc><c>Au</c><c> 1</c></rc><rc><c>O </c><c> 2</c></rc></set></array></atominfo>',
         '<calculation><scstep><energy><i name="e_fr_energy"> -1.0 </i></energy></scstep>',
         '<dos><i name="efermi"> %12.8f </i><total><array>' % EF,
         '<field>energy</field><field>total</field><field>integrated</field><set>']
    for s in range(ispin):
        X.append('<set comment="spin %d">' % (s + 1))
        X += ['<r> %10.4f %10.4f %10.4f </r>' % (e[i], tot[i, s], tin[i, s]) for i in range(NEDOS)]
        X.append('</set>')
    X.append('</set></array></total>')
    if not total_only:
        X.append('<partial><array><field>energy</field>' + ''.join('<field>%5s</field>' % o for o in ORBS) + '<set>')
        for a in range(2):
            X.append('<set comment="ion %d">' % (a + 1))
            for s in range(ispin):
                X.append('<set comment="spin %d">' % (s + 1))
                X += ['<r> %10.4f' % e[i] + ''.join(' %9.4f' % pd[a, s, o, i] for o in range(9)) + ' </r>'
                      for i in range(NEDOS)]
                X.append('</set>')
            X.append('</set>')
        X.append('</set></array></partial>')
    X += ['</dos></calculation>', '</modeling>']
    open(os.path.join(d, 'vasprun.xml'), 'w').write('\n'.join(X) + '\n')
    return pd


def _load(d, name):
    with open(os.path.join(d, name), 'rb') as f:
        return pickle.load(f)


@pytest.mark.parametrize('ispin', [1, 2])
def test_vasprun_fallback_matches_doscar(tmp_path, ispin):
    a, b = tmp_path / 'a', tmp_path / 'b'
    a.mkdir(); b.mkdir()
    pd = _write_pair(str(a), ispin)
    _write_pair(str(b), ispin)
    ra = common.get_dos(str(a) + '/')
    assert ra['source'] == 'DOSCAR' and not (a / 'dos_source.txt').exists()
    open(b / 'DOSCAR', 'w').write('')            # header-only / empty DOSCAR -> fallback
    rb = common.get_dos(str(b) + '/')
    assert rb['source'] == 'vasprun.xml' and (b / 'dos_source.txt').exists()
    da, db = np.array(_load(a, 'dos.pickle')), np.array(_load(b, 'dos.pickle'))
    assert da.shape == db.shape == (2, NEDOS - 1)
    np.testing.assert_allclose(da, db, atol=1e-9)
    np.testing.assert_allclose(da, pd.sum(axis=(1, 2))[:, 1:], atol=1e-9)
    np.testing.assert_allclose(_load(a, 'energy.pickle'), _load(b, 'energy.pickle'), atol=1e-9)
    if ispin == 2:
        sa, sb = _load(a, 'dos_spin.pickle'), _load(b, 'dos_spin.pickle')
        assert sa['source'] == 'DOSCAR' and sb['source'] == 'vasprun.xml'
        assert sa['orbitals'] == sb['orbitals'] and sa['symbols'] == sb['symbols'] == ['Au', 'O']
        for k in ('up', 'down', 'total_up', 'total_down', 'energy'):
            np.testing.assert_allclose(sa[k], sb[k], atol=1e-9)
    else:
        assert not (a / 'dos_spin.pickle').exists() and not (b / 'dos_spin.pickle').exists()


def test_vasprun_truncated_raises(tmp_path):
    _write_pair(str(tmp_path), 2)
    os.remove(tmp_path / 'DOSCAR')
    x = open(tmp_path / 'vasprun.xml').read()
    open(tmp_path / 'vasprun.xml', 'w').write(x[:len(x) // 2])
    with pytest.raises(ValueError, match='truncated'):
        common.get_dos(str(tmp_path) + '/')
    assert not (tmp_path / 'dos.pickle').exists()


def test_vasprun_total_only_raises(tmp_path):
    _write_pair(str(tmp_path), 2, total_only=True)
    os.remove(tmp_path / 'DOSCAR')
    with pytest.raises(ValueError, match='total DOS only'):
        common.get_dos(str(tmp_path) + '/')
