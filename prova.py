
# import logging



from dataclasses import dataclass
from dacite import from_dict, Config
from typing import Optional

@dataclass
class A:
    x: str
    y: int

class B(A):

    def __init__(self, x: str, y: int, z: Optional[int] = None):
        super().__init__(x, y)
        self.z = z

b = B(x='x', y=1, z=2)


print(b.x)
