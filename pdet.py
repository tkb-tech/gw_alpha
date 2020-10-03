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


m = [0.1,10.,20.,30.]
grid = len(m)-1
n = np.zeros((grid,grid))
v = np.zeros((grid,grid))
N = 100000
z_i = 0.00
z_f = 2.0
dz = z_f - z_i
#_alpha_ = 36.0/37.0 * (1.0 + tkb) * 0.9
_alpha_ = 36.0/37.0 * 0.8
H_0 = 70.0 #[(km/s)/Mpc]
c = 299792.0 #[km/s]
o_M = 0.3
o_k = 0.0
o_l = 0.7
s_s = 10000


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
    
    result = integrate.quad(g,0.0,1.0/(1.0+z))
    
    return result[0]/H_0


#describe about mass function
sigma = 0.6
m_c = 20.0 #unit is solar mass.
Const_in_R_th = 1.0

def f(m):
    return ( 1/ ( np.sqrt(2.0*np.pi) * sigma * m ) ) * math.exp( -( math.log(m/m_c) )**2.0 / ( 2.0 * sigma**2.0 ) )

"""
def f(m):
    if 5.0 < m and m < 50.0:
        return 1.0
    else:
        return 0.0
"""

def h(m):
    return m**(3.0/37.0) * f(m)


def R_th( m_1 ,m_2 ,z):
    return Const_in_R_th * t_c(z)**(-34.0/37.0) *( m_1 + m_2 )**( _alpha_ ) * h( m_1 ) * h( m_2 )

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
        return dVdz(z) * p( m_1 ,m_2 ,z ) * R_th( m_1 ,m_2 , z ) / ( 1.0 + z )

t1 = time.time()

for i in range(0,1):
    for j in range(0,1):
        """
        b= integrate.tplquad( lambda x,y,z: dR_det(x,y,z) ,z_i, z_f,lambda x:m[j], lambda x:m[j+1], lambda x,y:m[i], lambda x,y:m[i+1] )
        """
        
        sum_int = 0.0
        for l in range(0, N):
            sum_int += dR_det( random.uniform( m[i], m[i+1] ), random.uniform( m[j], m[j+1] ), random.uniform( z_i, z_f ) )
        
        b = (m[i+1]-m[i]) * (m[j+1]-m[j]) * dz * sum_int / N
        
        print b
        v[i][j] = b
        print 'Now (',i,',',j,')'

print v

t2 = time.time()
print 'the total time is'
print (t2-t1)/60.0
print 'minutes'

w = np.array([])

n1 = 10

n2 = 10

n3 = 100

sum_int_2 = 0.0

dh = (m[1] - m[0])/float(n1)
dz = (z_f - z_i)/float(n3)

t3 = time.time()

x = input()

if x == 1:
    for i in range(0,n1):
        print i
        print sum_int_2
        for j in range(0,n2):
            for k in range(0,n3):
                sum_int_2 = sum_int_2 + dz*dh*dh*(dR_det(m[0] + dh*(i+1),m[0] + dh*(j+1),z_i + dz*(k+1)))
elif x == 2:

    for i in range(0,n1):
        print i
        print sum_int_2
        for j in range(0,n2):
            for k in range(0,n3):
                sum_int_2 = sum_int_2 + dz*dh*dh*(dR_det(m[0] + dh*(i+0.5),m[0] + dh*(j+0.5),z_i + dz*(k+0.5)))

elif x == 3:
    for i in range(0,n1):
        print i
        print sum_int_2
        for j in range(0,n2):
            for k in range(0,n3):
                sum_int_2 = sum_int_2 + dz*dh*dh*(dR_det(m[0] + dh*(i),m[0] + dh*(j),z_i + dz*(k)))


            
t4 = time.time()


print sum_int_2

print ' '
print (t4-t3)/60.
print 'minutes'


