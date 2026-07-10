# common

Common utilities for constant-potential DFT with VASP/VASPsol++ and ASE,
shared across research groups (kpoint helpers, work-function/grand-canonical
energy parsing, and the `const_U_relax` / `const_U_dimer` / `const_U_FBL`
constant-potential drivers built on `set_pot`).

## Install

Either of these works; pick whichever fits how your group already runs jobs.

**Option A: PYTHONPATH** (no install step, matches existing usage)
```bash
export PYTHONPATH=/Users/jgauth32/PythonModules:$PYTHONPATH
```
Then `import common` from anywhere.

**Option B: pip install -e** (adds it to a venv's site-packages as an
editable package, useful if you don't want to manage PYTHONPATH)
```bash
pip install -e /Users/jgauth32/PythonModules/common
```

## Usage

```python
import common
from ase.calculators.vasp import Vasp
from ase.io import read

atoms = read('POSCAR')
calc = Vasp(xc='pbe', encut=400, kpts=(4, 4, 1), lsol=True)  # set your tags
common.const_U_relax(atoms, calc, desired_U=-0.5)  # V vs SHE
```

See `help(common)` or `help(common.const_U_relax)` for the full function
list and current keyword arguments (docstrings are the source of truth,
not this file).

## Examples

See `examples/` for runnable driver templates (`const_pot_relax.py`,
`const_pot_dimer.py`, `const_pot_fbl.py`) that call these functions with a
realistic `Vasp` calculator setup -- copy one into a job directory and edit
the potential, atom indices, and kpts for your system.
