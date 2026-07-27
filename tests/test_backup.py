"""Tests for handle_restart in common/__init__.py.

None of these touch VASP -- they exercise the backup_NN numbering, the
already-preserved check that lets the routine be called unconditionally,
and the restart staging for both relaxations and dimer runs, all against
fake VASP output files in a tmp_path.
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


def test_fresh_directory_is_left_alone(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _touch('POSCAR', 'starting geometry\n')

    assert common.handle_restart() is False
    assert not os.path.exists('backup_00')
    with open('POSCAR') as f:
        assert f.read() == 'starting geometry\n'


def test_backs_up_and_returns_path(rundir):
    path = common.handle_restart(dimer=True, quiet=True)

    assert path == 'backup_00'
    for name in ['INCAR', 'OUTCAR', 'vasprun.xml', 'CENTCAR', 'DIMCAR',
                 'MODECAR', 'NEWMODECAR', 'nelect_data.pkl']:
        assert os.path.isfile(os.path.join(path, name))
    # large files are deliberately excluded
    assert not os.path.exists(os.path.join(path, 'WAVECAR'))


def test_backup_preserves_contents(rundir):
    _touch('DIMCAR', 'step 1 curvature\n')
    path = common.handle_restart(dimer=True, quiet=True)

    with open(os.path.join(path, 'DIMCAR')) as f:
        assert f.read() == 'step 1 curvature\n'


def test_missing_files_are_skipped(rundir):
    os.remove('DIMCAR')
    os.remove('NEWMODECAR')

    path = common.handle_restart(dimer=True, quiet=True)

    assert not os.path.exists(os.path.join(path, 'DIMCAR'))
    assert os.path.isfile(os.path.join(path, 'OUTCAR'))


def test_relax_run_skips_the_dimer_files(rundir):
    path = common.handle_restart(quiet=True)

    assert not os.path.exists(os.path.join(path, 'DIMCAR'))
    assert os.path.isfile(os.path.join(path, 'CONTCAR'))


def test_extra_files_are_backed_up(rundir):
    _touch('iter03.traj')
    path = common.handle_restart(extra=['iter03.traj'], quiet=True)

    assert os.path.isfile(os.path.join(path, 'iter03.traj'))


def test_explicit_file_list_overrides_the_default(rundir):
    path = common.handle_restart(files=['CONTCAR', 'DIMCAR'], quiet=True)

    assert os.path.isfile(os.path.join(path, 'DIMCAR'))
    assert not os.path.exists(os.path.join(path, 'OUTCAR'))


def test_already_preserved_output_is_not_backed_up_twice(rundir):
    """A cleanly ended cycle already backed itself up -- do nothing."""
    assert common.handle_restart(dimer=True, quiet=True) == 'backup_00'
    assert common.handle_restart(dimer=True, quiet=True) is False
    assert not os.path.exists('backup_01')


def test_interrupted_run_is_preserved(rundir):
    """OUTCAR moved on after the last backup: the run was stopped."""
    common.handle_restart(dimer=True, quiet=True)
    _touch('OUTCAR', 'a longer OUTCAR from the run that got killed\n')

    assert common.handle_restart(dimer=True, quiet=True) == 'backup_01'


def test_output_never_preserved_counts_as_interrupted(rundir):
    """A directory from the older manual-backup workflow."""
    os.mkdir('backup_00')  # e.g. a partial hand-made backup

    assert common.handle_restart(dimer=True, quiet=True) == 'backup_01'


def test_a_match_in_any_backup_counts_not_only_the_newest(rundir):
    common.handle_restart(dimer=True, quiet=True)
    os.mkdir('backup_01')

    assert common.handle_restart(dimer=True, quiet=True) is False


def test_numbering_survives_a_gap(rundir):
    """Counting directories would collide with the existing backup_02."""
    os.mkdir('backup_00')
    os.mkdir('backup_02')

    assert common.handle_restart(quiet=True) == 'backup_03'


def test_unrelated_directories_are_ignored(rundir):
    os.mkdir('backup_old')
    os.mkdir('backups')

    assert common.handle_restart(quiet=True) == 'backup_00'


def test_prefix_is_configurable(rundir):
    assert common.handle_restart(prefix='cycle', quiet=True) == 'cycle_00'


def test_relax_run_continues_from_contcar(rundir):
    _touch('CONTCAR', 'relaxed geometry\n')
    _touch('MODECAR', 'old mode\n')
    _touch('NEWMODECAR', 'new mode\n')

    common.handle_restart(quiet=True)

    with open('POSCAR') as f:
        assert f.read() == 'relaxed geometry\n'
    # a relaxation has no mode to promote
    with open('MODECAR') as f:
        assert f.read() == 'old mode\n'


def test_dimer_run_continues_from_centcar_and_newmodecar(rundir):
    _touch('CENTCAR', 'dimer center\n')
    _touch('MODECAR', 'old mode\n')
    _touch('NEWMODECAR', 'new mode\n')

    path = common.handle_restart(dimer=True, quiet=True)

    with open('POSCAR') as f:
        assert f.read() == 'dimer center\n'
    with open('MODECAR') as f:
        assert f.read() == 'new mode\n'
    with open('oldMODECAR') as f:
        assert f.read() == 'old mode\n'
    # the backup keeps the mode the finished run *started* from
    with open(os.path.join(path, 'MODECAR')) as f:
        assert f.read() == 'old mode\n'


def test_mode_promotion_is_idempotent(rundir):
    _touch('MODECAR', 'mode\n')
    _touch('NEWMODECAR', 'mode\n')

    common.handle_restart(dimer=True, quiet=True)
    _touch('OUTCAR', 'output from the next cycle\n')
    common.handle_restart(dimer=True, quiet=True)

    with open('MODECAR') as f:
        assert f.read() == 'mode\n'


def test_wavecar_is_removed_by_default(rundir):
    common.handle_restart(dimer=True, quiet=True)

    assert not os.path.exists('WAVECAR')


def test_wavecar_is_kept_when_asked(rundir):
    """Between cycles of a running job the wavefunction is worth keeping."""
    common.handle_restart(dimer=True, rm_wavecar=False, quiet=True)

    assert os.path.isfile('WAVECAR')
