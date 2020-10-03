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
Const_in_R_th = 0.01*10**(1)

def f(m):
    return ( 1.0/ ( np.sqrt(2.0*np.pi) * sigma * m ) ) * math.exp( -( math.log(m/m_c) )**2.0 / ( 2.0 * sigma**2.0 ) )

def h(m):
    return m**(3.0/37.0) * f(m)


def R_th( m_1 ,m_2 ,z):
    if m_1<0.0:
        return 0.0
    elif m_2<0.0:
        return 0.0
    elif z < 0.0:
        return 0.0
    else:
        return Const_in_R_th * t_c(z)**(-34.0/37.0) *( m_1 + m_2 )**( 36.0 / 37.0 ) * h( m_1 ) * h( m_2 )


print 'input m_1'
m_1 = input()

print 'input m_2'
m_2 = input()

print 'input z'
z = input()

print R_th( m_1 ,m_2 ,z )

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

class Metropolis(object):
    
    def __init__(self, func, ndim, proposal_std=1., sample_rate=1):
        self.func = func
        self.ndim = ndim
        self.proposal_std = proposal_std
        self.sample_rate = sample_rate
    
    def __call__(self, sample_size):
        x = np.full(self.ndim ,50.0)
        samples = []
        for i in xrange(sample_size * self.sample_rate):
            x_new = np.random.normal(scale=self.proposal_std, size=self.ndim)
            x_new += x
            accept_prob = self.func(x_new) / self.func(x)
            if accept_prob > np.random.uniform():
                x = x_new
            if i % self.sample_rate == 0:
                samples.append(x)
        assert len(samples) == sample_size
        return np.array(samples)




def func(x):
    return R_th( x[0] ,x[1] ,x[2]/100.0 )

print "R_det sampling..."

s_s = 100

sampler = Metropolis(func, 3, proposal_std=15., sample_rate=15)
samples = sampler(sample_size=2*s_s)

print "mean\n", np.mean(samples, axis=0)
print "covariance\n", np.cov(samples, rowvar=False)

data_gw = np.array([])
data_m1 = np.array([])
data_m2 = np.array([])

for i in range (1,s_s):
    box_gw = samples[i-1]
    print box_gw
    
    data_m1 = np.append(data_m1 ,box_gw[0])
    
    data_m2 = np.append(data_m2 ,box_gw[1])
    
    data_gw = np.append(data_gw , box_gw)

plt.scatter(data_m1,data_m2)
plt.xlabel("m1")
plt.ylabel("m2")
plt.grid(True)
plt.show()

print samples



grid_1 = 0.00001
grid_2 = 0.00001
grid_3 = 0.00001
grid_4 = 0.00001

m = 30.0
dm = 30.0

for i in range(0,s_s - 1):
    data = samples[i+s_s]
    print data
    m1 = data[0]
    m2 = data[1]
    z = data[2]
    if m < m1 < m+dm and m < m2 < m+dm:
        grid_1 += 1
    elif m > m1 > m-dm and m < m2 < m+dm:
        grid_2 += 1
    elif m > m1 > m-dm and m > m2 > m-dm:
        grid_3 += 1
    elif m < m1 < m+dm and m > m2 > m-dm:
        grid_4 += 1



def deff_lnR(x1,x2,x3,x4):
    #m1 differential
    dh = dm
    dk = dm
    def f(a,b):
        return ( math.log( a ) - math.log( b ) )/dh
    def g(a,b,c,d):
        return ( f(a,b) - f(d,c) )/dk
    return g(x1,x2,x3,x4)


def alfa(a1,a2,a3,a4):
    return -((m+m)**2)*deff_lnR(a1,a2,a3,a4)

print grid_1,grid_2,grid_3,grid_4

print alfa(grid_1,grid_2,grid_3,grid_4)


