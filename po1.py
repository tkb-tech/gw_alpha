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
    y = o_M*(1.0+x)**3.0 - o_k*(1+x)**2.0 + o_l
    result = np.sqrt(y)
    return result


#describe comoving distance from redshhift x.
def D_c(x):
    
    def g(y):
        return 1.0/E(y)
    
    f = integrate.quad(g,0.0,x)
    
    return f[0]*c/H_0


def dVdz(z):
    
    x = (4.0*(np.pi)*c/H_0)*(D_c(z))*(D_c(z))/E(z)
    return x

#print 'please input redshift...'

#z = input()

#print dVdz(z)

#dV/dz is completed!


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
Const_in_R_th = 1

def f(m):
    return ( 1/ ( np.sqrt(2.0*np.pi) * sigma * m ) ) * math.exp( -( math.log(m/m_c) )**2.0 / ( 2.0 * sigma**2.0 ) )

def h(m):
    return m**(3.0/37.0) * f(m)


def R_th( m_1 ,m_2 ,z):
    return Const_in_R_th * t_c(z)**(-34.0/37.0) *( m_1 + m_2 )**( 36.0 / 37.0 ) * h( m_1 ) * h( m_2 )


#print 'input m_1'
#m_1 = input()

#print 'input m_2'
#m_2 = input()

#print 'input z'
#z = input()

#print R_th( m_1 ,m_2 ,z )

##R_th is completed!

#setting p_det

#print('please input m1')
#m1 = input()
#print('please input m2')
#m2 = input()
#print('please input z')
#z = input()

#p=gwdet.detectability()
#m1=10. # Component mass in Msun
#m2=10. # Component mass in Msun
#z=0.1  # Redshift

#print(p(m1,m2,z))  # Fraction of detectabile sources

#p_det is completed!

p = gwdet.detectability()

def dR_det( m_1 ,m_2 ,z ):
    if m_1 < 0.0:
        return 0.0
    elif m_2 < 0.0:
        return 0.0
    elif z < 0.0:
        return 0.0
    else:
        return dVdz(z) * p( m_1 ,m_2 ,z ) * R_th( m_1 ,m_2 , z ) / ( 1 + z )
#return 1 / ( 1 + z ) * 1 * p( m_1 ,m_2 ,z ) * R_th( m_1 ,m_2 , z )


def re_R_th(m1,m2,z):
    return dR_det( m1 ,m2 ,z )*( 1 + z )/dVdz( z )/(0.0001 + p( m1 ,m2 ,z ))




map_th = np.array([])

map_det = np.array([])

z = 0.3

for i in range(0,20):
    for j in range(0,20):
        map_th = np.append(map_th,R_th(0.001+5.0*i,0.001+5.0*j,z))
        print R_th(0.001+5.0*i,0.001+5.0*j,z)

for i in range(0,20):
    for j in range(0,20):
        map_det = np.append(map_det,dR_det(0.001+5.0*i,0.001+5.0*j,z))
        print re_R_th(0.001+5.0*i,0.001+5.0*j,z)



ax = plt.subplot(111)
im = ax.imshow(map_th.reshape(20,20))



# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)

plt.show()

ax = plt.subplot(111)
im = ax.imshow(map_det.reshape(20,20))



# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)

plt.show()

