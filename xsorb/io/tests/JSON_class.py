# from xsorb.adsorptiondata.adsorptionstructure import MoleculeRotation
from ase import Atoms


atoms = Atoms('H2O', positions=[[0, 0, 0], [0, 0.757, 0.587], [0, -0.757, 0.587]])

# aaa = MoleculeRotation(
#     atoms,
#     xrot="45.0",
#     yrot="30.0",
#     zrot="60",
#     mol_atom=1
# )

# # # Encode to JSON
# # json_data = encode(aaa)
# # print("Encoded JSON:", json_data)


# # # Decode from JSON
# # decodeddict = decode(json_data)
# # new_aaa = MoleculeRotation.fromdict(decodeddict)
# # print("Decoded object:", new_aaa)

# from ase.db import connect
# db = connect('test.json')

# db.write(atoms, data={'molecule_rotation': aaa})


# # Read the data back
# with connect('test.json') as db:
#     for row in db.select():
#         print("Row ID:", row.id)
#         xxx = row.data['molecule_rotation']
#         new_aaa = MoleculeRotation.fromdict(xxx)
#         print("Decoded object:", new_aaa)
#         print("Z Rotation:", new_aaa.zrot)

from dataclasses import dataclass, asdict

@dataclass
class Nested:
    """
    Example dataclass to demonstrate JSON serialization and deserialization.
    """
    nested_value: int

    def todict(self) -> dict:
        return self.__dict__

    @classmethod
    def fromdict(cls, data: dict) -> 'Nested':
        return cls(**data)

@dataclass
class prova:
    """
    Example dataclass to demonstrate JSON serialization and deserialization.
    """
    name: str
    value: int
    nested: Nested

    def todict(self) -> dict:
        return self.__dict__

    @classmethod
    def fromdict(cls, data: dict) -> 'prova':
        data['nested'] = Nested.fromdict(data['nested'])
        return cls(**data)


x = prova(
    name="example",
    value=42,
    nested=Nested(nested_value=100)
)

# print("Original Instance:", x)

# # Encode to dictionary
# encoded_dict = x.todict()
# print("Encoded Dictionary:", encoded_dict)

# # Decode from dictionary
# decoded_instance = prova.fromdict(encoded_dict)
# print("Decoded Instance:", decoded_instance)

from ase.db import connect
db = connect('test.json')

# db.write(atoms, bonds=[])

# for row in db.select():
#     print("Row ID:", row.id)
#     print("bonds:", row.bonds)

# db.write(atoms, data={'prova': x})


# # Read the data back
# for row in db.select():
#     print("Row ID:", row.id)
#     xxx = row.data['prova']
#     new_aaa = prova.fromdict(xxx)
#     print("Decoded object:", new_aaa)



from ase.constraints import FixCartesian

c = FixCartesian(a=0, mask=[True, False, True])

atoms.set_constraint(c)

db.write(atoms)