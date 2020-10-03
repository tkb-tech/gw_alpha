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
s_s = 100

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
Const_in_R_th = 0.01*10**(-10)

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

#class Metropolis(object):
    
#    def __init__(self, func, ndim, proposal_std=1., sample_rate=10):
#        self.func = func
#       self.ndim = ndim
#       self.proposal_std = proposal_std
#       self.sample_rate = sample_rate
    
#   def __call__(self, sample_size):
#        x = np.full(self.ndim ,30.0)
#        samples = []
#        for i in xrange(sample_size * self.sample_rate):
#            x_new = np.random.normal(scale=self.proposal_std, size=self.ndim)
#            x_new += x
#            accept_prob = self.func(x_new) / self.func(x)
#            if accept_prob > np.random.uniform():
#                x = x_new
#            if i % self.sample_rate == 0:
#               samples.append(x)
#       assert len(samples) == sample_size
#       return np.array(samples)




def func(x):
    return dR_det( x[0] ,x[1] ,x[2]/100.0 )

print "R_det sampling..."

fileobj = open( 'data_100000.text' , 'r' )

print fileobj

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


cho = []

for i in range(0,s_s):

    x = samlist[random.randrange(yososu/3)]
    cho.append(x)
samples = np.array(cho)
print samples

print "mean\n", np.mean(samples, axis=0)
print "covariance\n", np.cov(samples, rowvar=False)


data_gw = np.array([])
data_m1 = np.array([])
data_m2 = np.array([])

for i in range (1,99):
    box_gw = samples[random.randrange(s_s - 1)]
    
    data_m1 = np.append(data_m1 ,box_gw[0])
    
    data_m2 = np.append(data_m2 ,box_gw[1])
    
    data_gw = np.append(data_gw , box_gw)

plt.scatter(data_m1,data_m2)
plt.xlabel("m1")
plt.ylabel("m2")
plt.grid(True)
plt.show()

print samples

#samples1 = np.array([])

#for i in range (0,s_s):
#    samples1 = np.append(samples1,samples[(rate - 1)*s_s+i-1])


#print samples1

data = samples.T

kde = stats.gaussian_kde(data)

result_rate = np.array([])

for i in range(0,20):
    
    for j in range(0,20):
        
        result_rate = np.append(result_rate,10**4*kde([5.01*i+0.01,5.01*j+0.01,20])/dR_det(5.01*i+0.1,5.01*j+0.1,0.2))



ax = plt.subplot(111)
im = ax.imshow(result_rate.reshape(20,20))



# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)

plt.show()

#x, y = np.meshgrid(np.linspace(0, 100, 100), np.linspace(0, 100, 100))
#z = func(np.array([x, y]).reshape(2, -1).T).reshape(100, 100)
#plt.contour(x, y, z)
#plt.scatter(samples[:, 0], samples[:, 1])
#plt.xlim(-10, 10)
#plt.ylim(-10, 10)
#plt.show()


#N = 10*30*30*100
#sum = 0.0

#def poyo(m1,m2,z):
#    return dR_det( m1 ,m2 ,z ) * z

#for l in range(1, N):
#sum += dR_det( random.uniform( 1.0, 150.0 ) ,random.uniform( 1.0 , 150.0 ) ,random.uniform( 0.01,1.5 ))

#print ' denominator is ... '
#print 149 * 149 * 1.49 * sum / N


N = 30

dm_1 = 5

dm_2 = 5

#dz = 0.02

dz = 0.001

result = np.array( [] )

#for i in range(0,20):

#m_1 = 5.01 * ( i + 0.0001 )

#for j in range(0,20):

#m_2 = 5.01 * ( j + 0.0001 )


#z = 0.2

#sum = 0.0

#for l in range(1, N):
#sum += dR_det( random.uniform( m_1 - dm_1*0.5, m_1 + dm_1*0.5 ) ,random.uniform( m_2 - dm_2*0.5 ,m_2 + dm_2*0.5 ) ,random.uniform( z, z + dz ))

#print ' m1 = %f , m2 = %f ,z = %f ' % ( m_1 ,m_2 ,z)
#print dm_1 * dm_2 * dz * sum / N

#result = np.append( result, dm_1 * dm_2 * dz * sum / N )

#print result

#ax = plt.subplot(111)
#im = ax.imshow(result.reshape((20,20)))

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)

#plt.colorbar(im, cax=cax)

#plt.show()


def dR_z(m1,m2,z):
    return (1.0+z)*((10.0**4.0)*kde([m1,m2,100.0*z]))/(0.01+p(m1,m2,z))/dVdz(z)


def R_z(m1,m2):
    z_i = 0.3
    z_f = 0.6
    def po(z):
        return dR_z(m1,m2,z)
    answ = integrate.quad(po,z_i,z_f)
    return answ[0]

resul= np.array([])

for i in range(0,20):
    for j in range(0,20):
        resul = np.append(resul,dR_z(5.0*i+0.001,5.0*j+0.001,0.4))

ax = plt.subplot(111)
im = ax.imshow(resul.reshape((20,20)))

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)

plt.show()

def deff_lnR_z(m1,m2,z):
    #m1 differential
    dh = 15.0
    dk = 15.0
    if dR_z(m1 + dh,m2+dk,z)>0.0 and dR_z(m1 - dh,m2+dk,z)>0.0 and dR_z(m1+dh ,m2 - dk,z)>0.0 and dR_z(m1-dh ,m2-dk,z)>0.0:
        def f(m1,m2,z):
            return ( math.log( dR_z(m1 + dh,m2,z) ) - math.log( dR_z(m1 - dh,m2,z) ) )/(2*dh)
        def g(m1,m2,z):
            return (f(m1,m2 + dk,z)-f(m1,m2 - dk,z))/(2*dk)
        return g(m1,m2,z)
    else:
        return 0.0


def alfa_z(m1,m2,z):
    return -((m1+m2)**2)*deff_lnR_z(m1,m2,z)



def deff_lnR(m1,m2):
    #m1 differential
    dh = 5.0
    dk = 5.0
    def f(m1,m2):
        return ( math.log( R_z(m1 + dh,m2) ) - math.log( R_z(m1 - dh,m2) ) )/(2*dh)
    def g(m1,m2):
        return (f(m1,m2 + dk)-f(m1,m2 - dk))/(2*dk)
    return g(m1,m2)



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



def alfa(m1,m2):
    return -((m1+m2)**2)*deff_lnR(m1,m2)


def neo_alfa(m1,m2,z):
    dm = 5.0
    return math.log( dR_z(m1,m2,z) * dR_z(m1+dm,m2+dm,z) / dR_z(m1+dm,m2,z) / dR_z(m1,m2+dm,z) ) / math.log((m1+m2)*(m1+m2+2.0*dm)/(m1+m2+dm)**2.0)





for n in range(1000):
    print 'input z'
    z = input()
    ha = np.array([])
    for i in range(0,40):
        for j in range(0,40):
            wasd = 0.0
            m1 = 25.1 + i
            m2 = 25.1 + j
            wasd = alfa_z(m1,m2,z)
            ha = np.append(ha,wasd)


    ax = plt.subplot(111)
    im = ax.imshow(ha.reshape(40,40),cmap="inferno")



# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    plt.colorbar(im, cax=cax)

    plt.show()

