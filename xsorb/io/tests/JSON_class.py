from dataclasses import dataclass, asdict
from dacite import from_dict, Config
from typing import Union


def dict_without_none(data):
    '''Return a dictionary excluding keys with None values.'''
    return dict(x for x in data if x[1] is not None)


@dataclass
class NestedA:
    """
    Example dataclass to demonstrate JSON serialization and deserialization.
    """
    nested_value: int
    __objtype__ : str = 'NestedA'


@dataclass
class NestedB:
    """
    Example dataclass to demonstrate JSON serialization and deserialization.
    """
    nested_value: int
    possible_none: int | None = None
    __objtype__ : str = 'NestedB'


@dataclass
class Intermediate:
    '''
    Class that contains either NestedA or NestedB
    '''
    nested: NestedA | NestedB


@dataclass
class External:
    """
    External class that should be directly serializable.
    """
    name: str
    value: int
    intermediate: Intermediate


    def todict(self):
        return asdict(self, dict_factory=dict_without_none)

    @classmethod
    def fromdict(cls, data: dict):
        def union_type_hook(value):
            # Gestisce la union NestedA | NestedB
            if isinstance(value, dict) and '__objtype__' in value:
                objtype = value['__objtype__']
                if objtype == 'NestedA':
                    return from_dict(NestedA, value)
                elif objtype == 'NestedB':
                    return from_dict(NestedB, value)
            raise ValueError(f"Unknown objtype: {value.get('__objtype__', 'missing') if isinstance(value, dict) else 'not a dict'}")

        return from_dict(cls, data, config=Config(type_hooks={
            Union[NestedA, NestedB]: union_type_hook
        }))


x = External(
    name="example",
    value=42,
    intermediate=Intermediate(NestedB(nested_value=100))
)

print("Original Instance:", x)

# Encode to dictionary
encoded_dict = x.todict()
print("Encoded Dictionary:", encoded_dict)

# Decode from dictionary
decoded_instance = External.fromdict(encoded_dict)
print("Decoded Instance:", decoded_instance)


print("Equality check:", x == decoded_instance)



# from ase.db import connect
# from ase import Atoms
# atoms = Atoms('H2O', positions=[[0, 0, 0], [0, 0.757, 0.587], [0, -0.757, 0.587]])
# db = connect('test.json')

# db.write(atoms, data={'prova': x})


# # Read the data back
# for row in db.select():
#     print("Row ID:", row.id)
#     xxx = row.data['prova']
#     new_aaa = External.fromdict(xxx)
#     print("Decoded object:", new_aaa)