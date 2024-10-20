from ase import Atoms

atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 1]], cell=[3, 3, 3])

print(atoms.get_scaled_positions(wrap=False).flatten())