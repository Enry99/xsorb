from ase.data.pubchem import pubchem_atoms_conformer_search

from ase.io import write


conformers = pubchem_atoms_conformer_search(name='H2O')

print(f'Found {len(conformers)} conformers')
