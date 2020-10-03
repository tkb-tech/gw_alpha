import gwdet
import numpy as np
import scipy
from scipy import integrate
import math
import random

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

Const_in_R_th = 0.00000000001

def t_c(z):
    def f(a):
        y  = o_M/a - o_k + o_l*a*a
        result = np.sqrt(y)
        return result
    def g(a):
        result = 1.0/f(a)
        return result
    
    result = integrate.quad(g,0,1.0/(1+z))
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
    return 1 / ( 1 + z ) * dVdz(z) * p(m_1 ,m_2 ,z) * R_th( m_1 ,m_2 , z )


#print 'input m_1'
#m_1 = input()
#print 'input dm_1'
#dm_1 = input()

#print 'input m_2'
#m_2 = input()
#print 'input dm_2'
#dm_2 = input()

#print 'input z'
#z = input()
#print 'input dz'
#dz = input()
print 'input N'
N = input()

#print dVdz(z)

#print p(m_1 ,m_2 ,z)

#print R_th( m_1 ,m_2 , z )

#print dR_det( m_1 ,m_2 ,z )


dm_1 = 5

dm_2 = 5

dz = 0.02

for i in range(1,15):
    
    m_1 = dm_1 * ( i + 1 )
    
    for j in range(1,15):
        
        m_2 = dm_2 * ( j + 1 )
        
        for k in range(4,6):
            
            
            z = ( k + 1 )*dz
        
    
            sum = 0.0

            for l in range(1, N):
                sum += dR_det( random.uniform( m_1, m_1 + dm_1 ) ,random.uniform( m_2 ,m_2 + dm_2 ) ,random.uniform( z, z + dz ))


            print ' m1 = %f , m2 = %f ,z = %f ' % ( m_1 ,m_2 ,z)
            print dm_1 * dm_2 * dz * sum / N
