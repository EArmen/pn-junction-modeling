import numpy as np


def eps_max(num, ref):
    return float(np.max(np.abs(num - ref)))


def eps_rms(num, ref):
    return float(np.sqrt(np.mean((num - ref) ** 2)))