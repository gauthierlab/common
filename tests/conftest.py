"""Pytest configuration: make sure `import common` resolves to this repo.

/Users/jgauth32/PythonModules/common is the package directory itself
(__init__.py at the repo root), so the thing that needs to be on sys.path
is its *parent*, /Users/jgauth32/PythonModules -- exactly like the
PYTHONPATH setup documented in the README.
"""
import sys
import os

_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PYTHONMODULES = os.path.dirname(_PARENT)

if _PYTHONMODULES not in sys.path:
    sys.path.insert(0, _PYTHONMODULES)
