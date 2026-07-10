# TODO

## Split `__init__.py` into submodules (deferred, some day)

The package is a single ~1400-line `__init__.py` covering four distinct concerns.
Planned split:

- `_io.py` — grep helpers and VASP output parsers (`greplines`, `get_wf_*`, `get_chgcar`, ...)
- `_analysis.py` — workfunction/DOS/free-energy-diagram analysis (`get_dos`, `fed`, `get_omega`, ...)
- `_geometry.py` — cell/geometry/constraint utilities (`pos_swap`, `match_cell`, `reindex_atoms`, ...)
- `_constpot.py` — `set_pot` and the `const_U_*` drivers

Hard requirements for whoever does this:

1. **Zero caller breakage.** `__init__.py` must re-export every public function so
   `import common; common.some_func(...)` and `from common import greplines, get_wf_implicit`
   keep working unchanged for all groups.
2. **Keep the module-constant override pattern working.** `_she_U` / `_tolerance_U` stay
   defined in `__init__.py`, and functions must read them at call time (the
   `she_U=None` kwarg resolution already does this) so `common._she_U = 4.6`
   still takes effect after the split.
3. Run the test suite (`python3 -m pytest tests/`) before and after; it must stay green.
4. `get_interp` has a tab-indented body — move it verbatim, don't re-indent.

## Other deferred items

- Update the standalone sibling scripts one directory up (`const_pot_fbl.py`,
  `const_pot_dimer.py`) to be thin drivers like `examples/` — they still carry
  stale local copies of `set_pot`/`fmax`.
- Longer-term idea (separate project, Fortran): implement a native fixed-bond-length
  algorithm inside VTST so `const_U_FBL` could use in-VASP optimization and benefit
  from density/wavefunction extrapolation between steps. ASE-side constraints can't
  do this, and VASP's ICONST is MD-only (IBRION=0).
