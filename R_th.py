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

H_0 = 70.0 #[(km/s)/Mpc]
c = 299792.0 #[km/s]
o_M = 0.3
o_k = 0.0
o_l = 0.7

#create dV/dz(z) start!

def E(x):
    y = o_M*(1+x)**3 - o_k*(1+x)**2 + o_l
    result = np.sqrt(y)
    return result


#describe comoving distance from redshhift x.
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

#dV/dz is completed!


##describe R_th

#describe cosmic time at merger event.

Const_in_R_th = 0.01

def t_c(z):
    def f(a):
        y  = o_M/a - o_k + o_l*a*a
        result = np.sqrt(y)
        return result
    def g(a):
        result = 1.0/f(a)
        return result
    
    result = integrate.quad(g,0,1.0/(1.0+z))
    
    return result[0]/H_0



#print 'please write redshift...'
#z = input()
#print t_c(z)*3*(10**19)/(3.154*10**7)




#describe about mass function
sigma = 0.6
m_c = 20.0 #unit is solar mass.

def f(m):
    return ( 1/ ( np.sqrt(2.0*np.pi) * sigma * m ) ) * math.exp( -( math.log(m/m_c) )**2.0 / ( 2.0 * sigma**2.0 ) )

def h(m):
    return m**(3.0/37.0) * f(m)


def R_th( m_1 ,m_2 ,z):
    return Const_in_R_th * t_c(z)**(-34.0/37.0) *( m_1 + m_2 )**( 150.0 / 37.0 ) * h( m_1 ) * h( m_2 )


def deff_lnR_th(m1,m2,z):
        #m1 differential
    dh = 15.0
    dk = 15.0
    def f(m1,m2,z):
        return ( math.log( R_th(m1 + dh,m2,z) ) - math.log( R_th(m1 - dh,m2,z) ) )/(2*dh)
    def g(m1,m2,z):
        return (f(m1,m2 + dk,z)-f(m1,m2 - dk,z))/(2*dk)
    return g(m1,m2,z)
    
def alfa_th(m1,m2,z):
    return -((m1+m2)**2)*deff_lnR_th(m1,m2,z)

result = np.array([])
al = 0.0
z = 0.6
for i in range(0,20):
    m1 = 20.0 + i*3.0
    for j in range(0,20):
        m2 = 20.0 + j*3.0
        al = alfa_th(m1,m2,z)
        result = np.append(result,al)

ax = plt.subplot(111)
im = ax.imshow(result.reshape((20,20)))

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)

plt.show()

