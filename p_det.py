import gwdet
import numpy as np
import scipy
from scipy import integrate
import math
import random
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable



#p=gwdet.detectability()
#m1=10. # Component mass in Msun
#m2=10. # Component mass in Msun
#z=0.1  # Redshift

#print(p(m1,m2,z))  # Fraction of detectabile sources

#p_det is completed!

p = gwdet.detectability()

print p(31.27493985 ,34.56865202 ,94.41293787/100.0)

result = np.array( [] )

for i in range(0,100):
    for j in range(0,100):
        result = np.append(result ,p( 0.0 + i , 0.0 + j , 0.3 ))

ax = plt.subplot(111)
im = ax.imshow(result.reshape((100,100)),cmap="Greys")

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)

plt.show()
