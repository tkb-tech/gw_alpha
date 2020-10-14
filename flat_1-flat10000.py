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
s_s = 10000
m = [5.0, 18.0, 25.0, 30.0, 34.0, 38.0, 41.0, 44.0, 46.0, 48.0, 50.0]
#m = [0.0,15.0,20.0,23.0,27.0,32.0,38.0,44.0,53.0,65.0,100.0]
#m = [0.0, 15.0, 19.0, 25.0, 29.0, 33.0, 40.0, 46.0, 55.0, 65.0, 100.0]
i_a = 1
grid = len(m) - 1

def poisson(n,v):
    return (v**n) * (np.exp(-v)) / math.factorial(n)

def multi(n,v):
    result = math.factorial(total)
    for i in range(0,grid):
        for j in range(0,grid):
            result = result * (v[i][j]/total)**(n[i][j]) / math.factorial(n[i][j])
    return result

print 'R_det sampling...'


fileobj = open( 'alfa_1_100000_flat.text' , 'r' )
#fileobj = open( 'data_100000.text' , 'r' )
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

alfa_result = np.array([])
cho = []

for i in range(0,s_s):
    x = samlist[random.randrange(yososu/3)]
    cho.append(x)
    print[x[0],x[1],x[2]]
    print ','


samples = np.array(cho)




redshift_list = []

for i in range(0,s_s):
    if samples[i][0] > samples[i][1]:
        first = samples[i][0]
        second = samples[i][1]
        samples[i][0] = second
        samples[i][1] = first
    redshift_list.append(samples[i][2])

z_max = 0.01 * max(redshift_list)
z_min = 0.01 * min(redshift_list)

print("max redshhift is ",z_max,".")
print("min redshhift is ",z_min,".")

mean_redshift = np.mean(samples, axis=0)
print "mean\n", np.mean(samples, axis=0)
print "covariance\n", np.cov(samples, rowvar=False)

# samples are created!


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
sigma = 0.6
m_c = 20.0 #unit is solar mass.
Const_in_R_th = 1.0

"""
def f(m):
    return ( 1/ ( np.sqrt(2.0*np.pi) * sigma * m ) ) * math.exp( -( math.log(m/m_c) )**2.0 / ( 2.0 * sigma**2.0 ) )
"""

def f(m):
    if 5.0 < m and m < 50.0:
        return 1.0
    else:
        return 0.0


def h(m):
    return m**(3.0/37.0) * f(m)

res_kai = []

probability = []
probability_multinomial = []



for tkb in range(0,i_a):

    _alpha_ = (1.0 + tkb) * inpt

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
    N = 10000
    z_i = z_min
    z_f = z_max
    dz = z_f - z_i

    for i in range(0,grid):
        for j in range(0,grid):
            """
            b= integrate.tplquad( lambda x,y,z: dR_det(x,y,z) ,z_i, z_f,lambda x:m[j], lambda x:m[j+1], lambda x,y:m[i], lambda x,y:m[i+1] )
            """
            
            
            sum_int = 0.0
            for l in range(0, N):
                sum_int += dR_det( random.uniform( m[i], m[i+1] ), random.uniform( m[j], m[j+1] ), random.uniform( z_i, z_f ) )
            
            b = (m[i+1]-m[i]) * (m[j+1]-m[j]) * dz * sum_int / N
            
            
            """
            w = np.array([])

            n1 = 30

            n2 = 30

            n3 = 100

            sum_int_2 = 0.0

            dhi = (m[i+1] - m[i])/float(n1)
            dhj = (m[j+1] - m[j])/float(n2)
            dz = (z_f - z_i)/float(n3)

            for i1 in range(0,n1):
                for j1 in range(0,n2):
                    for k1 in range(0,n3):
                        #sum_int_2 = sum_int_2 + dz*dhi*dhj*(dR_det(m[i] + dhi*(i1+0.5),m[j] + dhj*(j1+0.5),z_i + dz*(k1+0.5)))
                        sum_int_2 = sum_int_2 + dz*dhi*dhj*(dR_det(m[i] + dhi*(i1),m[j] + dhj*(j1),z_i + dz*(k1)))
            
            b = sum_int_2
            """
            
            print b
            v[i][j] = b
            print 'Now (',i,',',j,')'

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
                if m[i] < m_1 and m_1 < m[i+1] and m[j] < m_2 and m_2 < m[j+1]:
                    n[i][j] = n[i][j] + 1.0
                    total = total + 1.0

    result = np.array([])
    al = 0.0

    for i in range(0,grid):
        for j in range(0,grid):
            al = v[i][j]
            result = np.append(result,al)


    ax = plt.subplot(111)
    im = ax.imshow(result.reshape((grid,grid)))
    
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    plt.colorbar(im, cax=cax)
    
    plt.show()


    kai_n = 0.0
    for i in range(0,grid):
        for j in range(i,grid):
            if i != j:
                
                v[i][j] = v[i][j] + v[j][i]
                n[i][j] = n[i][j] + n[j][i]
                v[j][i] = 0.0
            v[i][j] = v[i][j] * total / sum_v
            kai_n = kai_n + (n[i][j]-v[i][j])**2.0/v[i][j]
    """
    def kai(d):
        g = 0.0
        for i in range(0,grid):
            for j in range(0,grid):
                g = g + (d[i][j]-v[i][j])**2/v[i][j]
    """

    print 'number of total samples is'
    print total
    
    print 'the v[i][j]s are'
    print v

    print 'the n[i][j]s are'
    print n

    print 'when the alpha is '
    print _alpha_
    print 'the kai-value is...'
    print kai_n
    print ' '
    
    x = [_alpha_,kai_n]
    res_kai.append(x)
    """
    p_n_alpha = 1.0
    for i in range(0,grid):
        for j in range(0,grid):
            p_n_alpha = p_n_alpha * poisson(n[i][j],v[i][j])
    print 'the p_n_alpha by poisson is'
    print p_n_alpha

    baki = multi(n,v)
    print 'the p_n_alpha by multinomial is'
    print baki
    #probability.append([_alpha_,p_n_alpha])
    probability_multinomial.append([_alpha_,baki])
    print ' '
    """

t2 = time.time()

print res_kai

"""
print probability
print probability_multinomial
"""

print 'the mass function is ;( 1/ ( np.sqrt(2.0*np.pi) * sigma * m ) ) * math.exp( -( math.log(m/m_c) )**2.0 / ( 2.0 * sigma**2.0 ) ) .'
print 'masses that splits the plane are'
print m

print 'the total time is'
print (t2-t1)/60.0
print 'minutes'
