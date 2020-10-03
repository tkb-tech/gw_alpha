import gwdet
import numpy as np
import scipy
from scipy import integrate
from scipy import stats

from mpl_toolkits.mplot3d import Axes3D

import math
import random
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
s_s = 10000



print "R_det sampling..."

fileobj = open( 'da1a.text' , 'r' )

print fileobj

ld = np.array([])

for line in fileobj:
    val = float(line)
    #    print val
    ld = np.append(ld,val)
yososu = ld.shape[0]

fileobj.close()

print yososu

samlist = []

for i in range(0,yososu/3):
    x = np.array([])
    x = np.append(x,ld[3*i])
    x = np.append(x,ld[3*i+1])
    x = np.append(x,ld[3*i+2])
    samlist.append(x)

#print ' '
#print resu1t

cho = []

for i in range(0,s_s):
    x = samlist[i]
#    x = samlist[random.randrange(yososu/3)]
    cho.append(x)
samples = np.array(cho)
print samples

print "mean\n", np.mean(samples, axis=0)
print "covariance\n", np.cov(samples, rowvar=False)

