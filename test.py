from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor

atoms = read('CONTCAR')
aaa = AseAtomsAdaptor()
aaa.get_structure(atoms)
