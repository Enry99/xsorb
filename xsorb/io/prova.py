from ase import Atoms
import ase.db
with ase.db.connect('aaa.db') as db:
    atoms = Atoms('H2O', positions=[[0,0,0],[0,0,1.1],[0,0,1.1]])

    import numpy as np

    x = np.array([1,2,3])

    db.update(3, data={"y": {"y.2": 2}})


    for row in db.select():
        print(row.data)