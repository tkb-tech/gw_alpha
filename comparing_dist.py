print('please input the value of alpha...')
inpt = input()
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

#time measurement.
t1 = time.time()

H_0 = 70.0 #[(km/s)/Mpc]
c = 299792.0 #[km/s]
o_M = 0.3
o_k = 0.0
o_l = 0.7
#_alpha_ = 36.0/37.0
_alpha_ = 1.
s_s = 100
grid = 3
#m = [0.1,21.0,30.0,40.0,52.0,150.]
m = [1.0, 20.0, 35.0, 150.0]
#m = [1.0, 23., 44., 150.]
#m = [5.0,31.,41.,50.]
i_a = 1

def poisson(n,v):
    return (v**n) * (np.exp(-v)) / math.factorial(n)

def multi(n,v):
    result = math.factorial(total)
    for i in range(0,grid):
        for j in range(0,grid):
            result = result * (v[i][j]/total)**(n[i][j]) / math.factorial(n[i][j])
    return result


#create dV/dz(z) start!

def E(x):
    y = o_M*(1.0+x)**3.0 - o_k*(1+x)**2.0 + o_l
    result = np.sqrt(y)
    return result


#describe comoving distance from redshift x.
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

Const_in_R_th = 1.0


#func_name = 'log-normal'
func_name = 'flat'
function_name = []

def f(m):
    global function_name
    if func_name == 'log-normal':
        sigma = 0.6
        m_c = 20.0 #unit is solar mass.
        function_name = [['name','log_normal'],['sigma',sigma],['m_c',m_c]]
        return ( 1/ ( np.sqrt(2.0*np.pi) * sigma * m ) ) * math.exp( -( math.log(m/m_c) )**2.0 / ( 2.0 * sigma**2.0 ) )
    
    elif func_name == 'flat':
        m_max = 50.
        m_min = 5.
        function_name = [['name','flat'],['m_max',m_max],['m_min',m_min]]
        if m_min <= m and m <= m_max:
            return 1.0
        else:
            return 0.


def h(m):
    return m**(3.0/37.0) * f(m)

res_kai = []

probability = []
probability_multinomial = []

p = gwdet.detectability()

for tkb in range(0,i_a):


    def R_th( m_1 ,m_2 ,z):
        return Const_in_R_th * t_c(z)**(-34.0/37.0) *( m_1 + m_2 )**( _alpha_ ) * h( m_1 ) * h( m_2 )

    def dR_det_p_1( m_1 ,m_2 ):
        if m_1 <= 0.0:
            return 0.0
        elif m_2 <= 0.0:
            return 0.0
        else:
            return ( m_1 + m_2 )**( _alpha_ ) * h( m_1 ) * h( m_2 )
    
    def dR_det( m_1 ,m_2 ,z ):
        if m_1 <= 0.0:
            return 0.0
        elif m_2 <= 0.0:
            return 0.0
        elif z <= 0.0:
            return 0.0
        else:
            return dVdz(z) * p( m_1 ,m_2 ,z ) * R_th( m_1 ,m_2 , z ) / ( 1.0 + z )
            
            
    _alpha_ = (1.0 + tkb) * inpt
    result = np.array([])
    N_int = 100
    zi = 0.01
    zf =2.0
    sum_res = 0.
    for i in range(0,100):
        m1 = i * 1.
        for j in range(0,100):
            m2 = j * 1.
            sum_integral = 0.
            for k in range(0,N_int):
                z = zi + (zf - zi) / N_int * k
                sum_integral = sum_integral + (zf - zi)/N_int * dR_det(m1,m2,z)
            
            result = np.append( result, sum_integral )
            print '(',i,',',j,')'
            sum_res = sum_res + sum_integral
            
    for i in range(0,10000):
        result[i] = result[i]/sum_res
    
    
        
        
    
    ax = plt.subplot(111)
    im = ax.imshow(result.reshape((100,100)))
    
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    plt.colorbar(im, cax=cax)
    
    plt.show()
    
    
    
    

    v = np.zeros((grid,grid))
    N = 100000
    #    z_i = z_min
    #    z_f = z_max
    #    dz = z_f - z_i
    
    sum = 0.
    result = np.array([])
    for i in range(0,grid):
        for j in range(0,grid):
            b= integrate.dblquad( lambda x,y: dR_det_p_1(x,y) ,m[j], m[j+1], lambda x:m[i], lambda x:m[i+1] )
            
            print b
            v[i][j] = b[0]
            print 'Now (',i,',',j,')'
            sum = sum + v[i][j]

    for i in range(0,grid):
        for j in range(0,grid):
            result = np.append( result, v[i][j]/sum )

    ax = plt.subplot(111)
    im = ax.imshow(result.reshape((grid,grid)))
    
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    plt.colorbar(im, cax=cax)
    
    plt.show()
    
    print 'the v[i][j]s are'
    print result.reshape((grid,grid))

    print 'when the alpha is '
    print _alpha_
    
    
    print('please input the value of alpha...')
    inpt = input()
    _alpha_ = (1.0 + tkb) * inpt

    v = np.zeros((grid,grid))
    N = 100000
    
    z_i = 0.1
    z_f = 2.0
    dz = z_f - z_i
    
    sum = 0.
    result = np.array([])
    for i in range(0,grid):
        for j in range(0,grid):
            #b= integrate.dblquad( lambda x,y: dR_det(x,y) ,m[j], m[j+1], lambda x:m[i], lambda x:m[i+1] )
            
            
            sum_int = 0.0
            for l in range(0, N):
                sum_int += dR_det( random.uniform( m[i], m[i+1] ), random.uniform( m[j], m[j+1] ), random.uniform( z_i, z_f ) )
            
            b = (m[i+1]-m[i]) * (m[j+1]-m[j]) * dz * sum_int / N
            
            print b
            v[i][j] = b
            print 'Now (',i,',',j,')'
            sum = sum + v[i][j]

    for i in range(0,grid):
        for j in range(0,grid):
            result = np.append( result, v[i][j]/sum )

    ax = plt.subplot(111)
    im = ax.imshow(result.reshape((grid,grid)))
    
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    plt.colorbar(im, cax=cax)
    
    plt.show()
    
    print 'the v[i][j]s are'
    print result.reshape((grid,grid))

    print 'when the alpha is '
    print _alpha_



t2 = time.time()

print res_kai

"""
print probability
print probability_multinomial
"""

print 'the mass function is',function_name
print 'masses that splits the plane are'
print m

print 'the total time is'
print (t2-t1)/60.0
print 'minutes'
