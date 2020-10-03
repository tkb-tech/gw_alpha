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
"""
#time measurement.
t1 = time.time()
"""

H_0 = 70.0 #[(km/s)/Mpc]
c = 299792.0 #[km/s]
o_M = 0.3
o_k = 0.0
o_l = 0.7
_alpha_ = 36.0/37.0
s_s = 121
m = [0.0,10.0,20.0,30.0,100.0]
grid = len(m) - 1
i_a = 1

def poisson(n,v):
    return (v**n) * (np.exp(-v)) / math.factorial(n)

def multi(n,v):
    result = math.factorial(total)
    for i in range(0,grid):
        for j in range(i,grid):
            result = result * (v[i][j]/total)**(n[i][j]) / math.factorial(n[i][j])
    return result


print "R_det sampling..."


"""
fileobj = open( 'data_100000.text' , 'r' )
#fileobj = open( 'alfa_4_1000.text' , 'r' )
#fileobj = open( 'alfa_4_100000.text' , 'r' )

#print fileobj

ld = np.array([])

for line in fileobj:
    val = float(line)
    #    print val
    ld = np.append(ld,val)

yososu = ld.shape[0]

fileobj.close()

print yososu


samlist = []

for i in range(0,yososu/3):
    x = np.array([])
    x = np.append(x,ld[3*i])
    x = np.append(x,ld[3*i+1])
    x = np.append(x,ld[3*i+2])
    samlist.append(x)

#print ' '
#print resu1t

alfa_result = np.array([])
cho = []

for i in range(0,s_s):
    x = samlist[random.randrange(yososu/3)]
    cho.append(x)
samples = np.array(cho)
"""
samples = [
     [ 23.2 , 2.59  ,5.1],
     [ 29.7 , 8.4  ,15.0],
     [ 39.6 , 29.4 ,34.0],
     [ 35.5 , 26.8 ,20.0],
     [ 30.5 , 25.3 ,11.0],
     [ 35.0 , 23.8 ,20.0],
     [ 50.2 , 34.0 ,49.0],
     [ 11.0 , 7.6 , 7.00],
     [ 30.8 , 20.0 ,20.0],
     [ 13.7 , 7.7 , 9.00],
     [ 23.2 ,13.6,  21.0],
     [ 35.6 ,30.6,   9.0]
           ]


for i in range(0,s_s):
    if samples[i][0] > samples[i][1]:
        first = samples[i][0]
        second = samples[i][1]
        samples[i][0] = second
        samples[i][1] = first

print samples
mean_redshift = np.mean(samples, axis=0)
print "mean\n", mean_redshift
print "covariance\n", np.cov(samples, rowvar=False)



#print samples[5][1]

# samples are created!



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
Const_in_R_th = 1.0**(-10)


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

"""
def h(m):
    return m**(-1.0/7.0) * f(m)
"""
res_kai = []

probability = []
probability_multinomial = []


t1 = time.time()

for tkb in range(0,i_a):

    _alpha_ = 36.0/37.0 * 1

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

    n = np.zeros((grid,grid))
    v = np.zeros((grid,grid))
    N = 100000
    z_i = 0.05
    z_f = 0.5
    dz = z_f - z_i

    for i in range(0,grid):
        for j in range(0,grid):
            sum_int = 0.0
            """
            b= integrate.tplquad( lambda x,y,z: dR_det(x,y,z) ,z_i, z_f,lambda x:m[j], lambda x:m[j+1], lambda x,y:m[i], lambda x,y:m[i+1] )
            """
            
            for l in range(0, N):
                sum_int += dR_det( random.uniform( m[i], m[i+1] ), random.uniform( m[j], m[j+1] ), random.uniform( z_i, z_f ) )
            
            b = (m[i+1]-m[i]) * (m[j+1]-m[j]) * dz * sum_int / N
            
            print b
            v[i][j] = b


    n = np.zeros((grid,grid))
    total = 0.0
    sum_v = 0.0

    for i in range(0,grid):
        for j in range(0,grid):
            sum_v = sum_v + v[i][j]
            for k in range(0,s_s):
                m_1 = samples[k][0]
                m_2 = samples[k][1]
                z = samples[k][2]
                if m[i] < m_1 and m_1 <= m[i+1] and m[j] < m_2 and m_2 <= m[j+1]:
                    n[i][j] = n[i][j] + 1.0
                    total = total + 1.0

    result = np.array([])
    al = 0.0

    for i in range(0,grid):
        for j in range(0,grid):
            al = v[i][j]
            result = np.append(result,al)

    """
    ax = plt.subplot(111)
    im = ax.imshow(result.reshape((grid,grid)))
    
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    plt.colorbar(im, cax=cax)
    
    plt.show()
    """

    def kai(d):
        g = 0.0
        for i in range(0,grid):
            for j in range(i,grid):
                g = g + (d[i][j]-v[i][j])**2/v[i][j]
        return g


    for i in range(0,grid):
        for j in range(i,grid):
            if i != j:
                v[i][j] = v[i][j] + v[j][i]
                n[i][j] = n[i][j] + n[j][i]
                v[j][i] = 0.0
            v[i][j] = v[i][j] * total / sum_v

    chi_value_0 = kai(n)
    print 'the chi-value is'
    print chi_value_0
    print ' '

    print 'number of total samples is'
    print total

    print 'the v[i][j]s are'
    print v

    print 'the n[i][j]s are'
    print n

    print 'when the alpha is '
    print _alpha_

    print ' '

    p_value = 0.0
    p_total = 0.0
    cnt = 0

    def cal_p(N,chi,i,j,k):
        if i == grid - 1 and j == grid - 1:
            k[i][j] = N
            chi = kai(k)
            y = multi(k,v)
            global p_total
            global cnt
            cnt = cnt + 1
            p_total = p_total + y
            #print k
            #print chi
            
            if chi > chi_value_0:
                global p_value
                p_value = p_value + y

        else:
            for q in range(0,N+1):
                k[i][j] = q
                if j < grid - 1:
                    cal_p(N-q,chi,i,j+1,k)
                else:
                    cal_p(N-q,chi,i+1,i+1,k)

    baki = multi(n,v)
    print 'the p_n_alpha by multinomial is'
    print baki
    print ' '

    print ' '
    print 'Now calculating the p-value for this samples.'
    N = 12
    chi = 0.0
    k = np.zeros((grid,grid))
    p_value = 0.0
    p_total = 0.0
    cal_p(N,chi,0,0,k)
    print 'the p-value is'
    print p_value
    print 'the normalization factor is'
    print p_total
    print 'the number of cases is'
    print cnt
    probability_multinomial.append([_alpha_,baki,p_value,p_total])

t2 = time.time()
print probability_multinomial
print 'the mass function is ;( 1/ ( np.sqrt(2.0*np.pi) * sigma * m ) ) * math.exp( -( math.log(m/m_c) )**2.0 / ( 2.0 * sigma**2.0 ) ) .'
print 'masses that splits the plane are'
print m

print 'the total time is'
print (t2-t1)/60.0
print 'minutes'
