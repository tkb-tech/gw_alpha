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
import time


grid = 100
z = 0.1


p = gwdet.detectability(psd = 'aLIGOEarlyHighSensitivityP1200087')
#p = gwdet.detectability()

result = np.array([])
for i in range(0,grid):
    for j in range(0,grid):
        result = np.append(result,p(i,j,z))

ax = plt.subplot(111)
im = ax.imshow(result.reshape((grid,grid)))

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)

plt.show()
