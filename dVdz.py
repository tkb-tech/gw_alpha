import gwdet
import numpy as np
import scipy
from scipy import integrate
import math
import matplotlib.pyplot as plt


H_0 = 70.0 #[(km/s)/Mpc]
c = 299792.0 #[km/s]
o_M = 0.3
o_k = 0.0
o_l = 0.7


def E(x):
    y = o_M*(1+x)**3 - o_k*(1+x)**2 + o_l
    result = np.sqrt(y)
    return result

def D_c(x):

    def g(y):
        return 1/E(y)
    
    f = integrate.quad(g,0,x)
    
    return f[0]*c/H_0


def dVdz(z):

    x = (4*(np.pi)*c/H_0)*(D_c(z))*(D_c(z))/E(z)
    return x

#print 'please input redshift...'

#z = input()

#print dVdz(z)

z = np.arange(0 ,1 ,0.01)

y = np.array( [] )

for i in range (0,100):
    y = np.append( y ,dVdz(z[i]) )

plt.plot(z,y)
plt.show()



