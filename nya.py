import numpy as np

def func(x):
    return np.exp(-0.5 * x[0]**2.0 * x[1]**2.0 / 5.)

v = np.array([1.0,5.0])

print func(v)
