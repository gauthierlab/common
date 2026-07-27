"""Tests for the backup helpers in common/__init__.py.

None of these touch VASP -- they exercise the backup_NN numbering, the
skip-missing-files behavior, and handle_restart's dimer restart shuffle,
all against fake VASP output files in a tmp_path.
"""
import os

import pytest

import common


def _touch(name, text=None):
    with open(name, 'w') as f:
        f.write(text if text is not None else name + '\n')


@pytest.fixture
def rundir(tmp_path, monkeypatch):
    """A tmp cwd holding a minimal set of fake dimer outputs."""
    monkeypatch.chdir(tmp_path)
    for name in ['INCAR', 'KPOINTS', 'POSCAR', 'CONTCAR', 'OUTCAR', 'OSZICAR',
                 'XDATCAR', 'vasprun.xml', 'opt.log', 'nelect_data.pkl',
                 'CENTCAR', 'DIMCAR', 'MODECAR', 'NEWMODECAR', 'WAVECAR']:
        _touch(name)
    return tmp_path


def test_make_backup_copies_files_and_returns_path(rundir):
    path = common.make_backup(common._DIMER_BACKUP_FILES, quiet=True)

    assert path == 'backup_00'
    assert os.path.isdir(path)
    for name in common._DIMER_BACKUP_FILES:
        assert os.path.isfile(os.path.join(path, name))
    # large files are deliberately excluded
    assert not os.path.exists(os.path.join(path, 'WAVECAR'))


def test_make_backup_preserves_contents(rundir):
    _touch('DIMCAR', 'step 1 curvature\n')
    path = common.make_backup(common._DIMER_BACKUP_FILES, quiet=True)

    with open(os.path.join(path, 'DIMCAR')) as f:
        assert f.read() == 'step 1 curvature\n'


def test_make_backup_skips_missing_files(rundir):
    os.remove('DIMCAR')
    os.remove('NEWMODECAR')

    path = common.make_backup(common._DIMER_BACKUP_FILES, quiet=True)

    assert not os.path.exists(os.path.join(path, 'DIMCAR'))
    assert os.path.isfile(os.path.join(path, 'OUTCAR'))


def test_make_backup_extra_files(rundir):
    _touch('iter03.traj')
    path = common.make_backup(common._BACKUP_FILES, extra=['iter03.traj'],
                              quiet=True)

    assert os.path.isfile(os.path.join(path, 'iter03.traj'))


def test_make_backup_increments(rundir):
    paths = [common.make_backup(quiet=True) for _ in range(3)]
    assert paths == ['backup_00', 'backup_01', 'backup_02']


def test_make_backup_numbering_survives_a_gap(rundir):
    """Counting directories (the old behavior) would collide with backup_02."""
    os.mkdir('backup_00')
    os.mkdir('backup_02')

    assert common.make_backup(quiet=True) == 'backup_03'


def test_make_backup_ignores_unrelated_directories(rundir):
    os.mkdir('backup_old')
    os.mkdir('backups')

    assert common.make_backup(quiet=True) == 'backup_00'


def test_make_backup_explicit_path(rundir):
    path = common.make_backup(path='snapshot', quiet=True)

    assert path == 'snapshot'
    assert os.path.isfile('snapshot/OUTCAR')


def test_make_backup_refuses_to_overwrite(rundir):
    os.mkdir('snapshot')
    with pytest.raises(OSError):
        common.make_backup(path='snapshot', quiet=True)


def test_handle_restart_returns_false_without_outcar(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert common.handle_restart() is False


def test_handle_restart_relax_restarts_from_contcar(rundir):
    _touch('CONTCAR', 'relaxed geometry\n')

    path = common.handle_restart()

    assert os.path.isfile(os.path.join(path, 'OUTCAR'))
    with open('POSCAR') as f:
        assert f.read() == 'relaxed geometry\n'
    assert not os.path.exists('WAVECAR')
    # dimer files are not part of the plain-relaxation backup
    assert not os.path.exists(os.path.join(path, 'DIMCAR'))


def test_handle_restart_dimer_restarts_from_centcar(rundir):
    _touch('CENTCAR', 'dimer center\n')
    _touch('MODECAR', 'old mode\n')
    _touch('NEWMODECAR', 'new mode\n')

    path = common.handle_restart(dimer=True)

    assert os.path.isfile(os.path.join(path, 'DIMCAR'))
    with open('POSCAR') as f:
        assert f.read() == 'dimer center\n'
    with open('MODECAR') as f:
        assert f.read() == 'new mode\n'
    with open('oldMODECAR') as f:
        assert f.read() == 'old mode\n'
    # the backup keeps the mode the finished run *started* from
    with open(os.path.join(path, 'MODECAR')) as f:
        assert f.read() == 'old mode\n'
