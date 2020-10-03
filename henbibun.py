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

##describe R_th

#describe cosmic time at merger event.



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
Const_in_R_th = 0.01*10**(-10)

def f(m):
    return ( 1/ ( np.sqrt(2.0*np.pi) * sigma * m ) ) * math.exp( -( math.log(m/m_c) )**2.0 / ( 2.0 * sigma**2.0 ) )

def h(m):
    return m**(3.0/37.0) * f(m)


def R_th( m_1 ,m_2 ,z):
    return Const_in_R_th * t_c(z)**(-34.0/37.0) *( m_1 + m_2 )**( 36.0 / 37.0 ) * h( m_1 ) * h( m_2 )



result = np.array([])

for i in range(0,20):
    for j in range(0,20):
        result = np.append(result,R_th(5.0*i+0.001,5.0*j+0.001,0.2))

ax = plt.subplot(111)
im = ax.imshow(result.reshape((20,20)))

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)

plt.show()


def R_th_z(m1,m2):
    def q(z):
        return R_th(m1,m2,z)
    ans = integrate.quad(q,0.1,1.0)
    return ans[0]

def deff_lnR(m1,m2):
    #m1 differential
    dh = 20.0
    dk = 20.0
    def f(m1,m2):
        return ( math.log( R_th_z(m1 + dh,m2) ) - math.log( R_th_z(m1 - dh,m2) ) )/(2*dh)
    def g(m1,m2):
        return (f(m1,m2 + dk)-f(m1,m2 - dk))/(2*dk)
    return g(m1,m2)

def neo_alfa(m1,m2,z):
    dm = 5.0
    return math.log( R_th(m1,m2,z) * R_th(m1+dm,m2+dm,z) / R_th(m1+dm,m2,z) / R_th(m1,m2+dm,z) ) / math.log((m1+m2)*(m1+m2+2*dm)/(m1+m2+dm)**2)

def alfa(m1,m2):
    return -((m1+m2)**2)*deff_lnR(m1,m2)

print 'input m1'
m1 = input()

print 'input m2'
m2 = input()

z = 0.3

print R_th(m1,m2,0.2)
print R_th(m1+0.01,m2,0.2)
print R_th_z(m1,m2)
print R_th_z(m1+0.01,m2)
print neo_alfa(m1,m2,z)

print alfa(m1,m2)








