import numpy as np
import scipy
from scipy import integrate
from scipy import stats

from mpl_toolkits.mplot3d import Axes3D

import math
import random
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

n = 10000
grid = 55

"""
p1 = 1./9.
p2 = 1./9.
p3 = 1./27.
p4 = 1./27.
p5 = 1./27.
p6 = 1./15.
p7 = 19./150.
p8 = 1./150.
p9 = 1./15.
p10 =1./15.
p11 =1./15.
p12 =1./15.
p13 =1./15.
p14 =1./15.
p15 =1./15.

pro = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15]
"""

pro = np.array([])

for i in range(0,grid):
    x = 1./grid
    pro = np.append(pro,x)


def kai(x):
    result = 0.0
    for i in range(0,len(pro)):
        result = result + (x[i]-n*pro[i])**2/(n*pro[i])
    return result

def sum(x):
    hoge = 0.
    for i in range(0,len(x)):
        hoge = hoge + x[i]
    return hoge



print (pro)

s_s = 100000

x = np.random.multinomial(n,pro,size = s_s)

print (x)


k = np.array([])

for i in range(0,s_s):
    k = np.append(k,kai(x[i]))

print (k)

y = np.random.chisquare(len(pro)-1,s_s)

print (y)

print (sum(pro))

plt.hist(k, bins=50,range=(0,30.0),color='red',alpha=0.7)
plt.hist(y, bins=50,range=(0,30.0),color='blue',histtype='step')
plt.show()





