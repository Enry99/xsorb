from chgnet.model import CHGNet
from chgnet.model import CHGNetCalculator


def NNloader():
    model = CHGNet.load()
    calculator = CHGNetCalculator(model)
    return calculator
