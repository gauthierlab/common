# examples

These are templates, not turnkey jobs -- copy the one you need into an
actual job directory and edit the potential (`pot_des`), any atom indices
(`o_ind`/`h_ind`/`ind1`/`ind2`), `kpts`, and other VASP tags for your system
before submitting.

- `const_pot_relax.py` -- constant-potential geometry relaxation
  (`common.const_U_relax`)
- `const_pot_dimer.py` -- constant-potential Dimer transition-state search
  (`common.const_U_dimer`)
- `const_pot_fbl.py` -- constant-potential fixed-bond-length optimization
  (`common.const_U_FBL`)

They exist so groups stop re-deriving `set_pot`/`fmax`/checkpointing logic
locally in every job directory -- the actual optimization logic lives once,
in `common/__init__.py`.
