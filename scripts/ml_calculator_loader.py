# This is an example of how to load an ase.calculator for a machine learning model
# (here CHGNet), but can be many others, such as MACE, or DeepMD, etc. (basically
# any potential that implements the ase.calculator interface. You can even use the
# Lennard-Jones potential if you want :) )

# The only requirement is that this script defines the function NNloader() that
# returns the calculator object. This function will be called by the main script.

from chgnet.model import CHGNet
from chgnet.model import CHGNetCalculator


def NNloader():
    model = CHGNet.load()
    calculator = CHGNetCalculator(model)
    return calculator
