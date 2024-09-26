from pathlib import Path

p = Path('LICENSE')

#print(p.is_file())

import numpy as np


x = np.array([{"a": 1, "b": 2}, {"a": 3, "b": 4}])

np.save("test.npy", x)

#y = np.load("test.npy")
#print(y)