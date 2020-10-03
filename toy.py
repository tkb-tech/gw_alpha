import numpy as np
import scipy
from scipy import integrate
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D
import math
import random
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable



grid = 3

m = [5.0,32.0,49.0,100]


def dR_det(x,y,z):
    return x**2+y+z**3
    
n = np.zeros((grid,grid))
v = np.zeros((grid,grid))
N = 100000
z_i = 0.05
z_f = 1.5
dz = z_f - z_i

b= integrate.tplquad( lambda x,y,z: dR_det(x,y,z) ,0.0, 0.01,lambda x:-100.0, lambda x:100.0, lambda x,y:0.0, lambda x,y:100.0 )
print b
