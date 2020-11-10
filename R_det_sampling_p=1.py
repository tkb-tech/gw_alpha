import gwdet
import numpy as np
import scipy
from scipy import integrate
import math
import random
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

H_0 = 70.0 #[(km/s)/Mpc]
c = 299792.0 #[km/s]
o_M = 0.3
o_k = 0.0
o_l = 0.7
s_s = 100000

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



"""
m_max = 50.
m_min = 5.
def f(m):
    if m_min < m and m < m_max:
        return 1.0
    else:
        return 0.0
"""

def h(m):
    return m**(3.0/37.0) * f(m)


def R_th( m_1 ,m_2 ,z):
    return Const_in_R_th * t_c(z)**(-34.0/37.0) * ( m_1 + m_2 )**( 36.0 / 37.0 ) * h( m_1 ) * h( m_2 )



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

#p = gwdet.detectability()

"""
p = 1.0

def dR_det( m_1 ,m_2 ,z ):
    if m_1 < 0.0:
        return 0.0
    elif m_2 < 0.0:
        return 0.0
    elif z < 0.0:
        return 0.0
    else:
        return 1 / ( 1 + z ) * dVdz(z) * p * R_th( m_1 ,m_2 , z )



p_det = gwdet.detectability()

def dR_det_pdet( m_1 ,m_2 ,z ):
    if m_1 < 0.0:
        return 0.0
    elif m_2 < 0.0:
        return 0.0
    elif z < 0.0:
        return 0.0
    else:
        return 1 / ( 1 + z ) * dVdz(z) * p_det( m_1 ,m_2 ,z ) * R_th( m_1 ,m_2 , z )




R_list = np.array([])
z_list = np.array([])

for i in range(0,100):
    z = i * 1.
    R_list = np.append( R_list, dR_det(20.,20.,z) )
    z_list = np.append( z_list, z )
    

R_list_p = np.array([])
z_list_p = np.array([])

for i in range(0,100):
    z = i * 0.01
    R_list_p = np.append( R_list_p, dR_det_pdet(20.,20.,z) )
    z_list_p = np.append( z_list_p, z )


plt.scatter(z_list,R_list,s = 0.01)
plt.plot(z_list,R_list)
plt.scatter(z_list_p,R_list_p,s = 0.01)
plt.plot(z_list_p,R_list_p,c = 'red')
plt.ylabel("R_det,p = 1")
plt.xlabel("z")
plt.yscale('log')
plt.grid(True)
plt.show()

plt.cla()
plt.clf()
"""

def dR_det(m_1,m_2):
    if m_1 <= 0.0:
        return 0.0
    elif m_2 <= 0.0:
        return 0.0
    else:
        return ( m_1 + m_2 )**( 1. ) * h( m_1 ) * h( m_2 )



#calculate the average of m1,m2 in dR_det.

integ1 = integrate.dblquad( lambda x,y: x*dR_det(x,y) , 3. , 60. , lambda x: 3 , lambda x:60. )

integ2 = integrate.dblquad( lambda x,y: y*dR_det(x,y) , 3. , 60. , lambda x: 3 , lambda x:60. )

denomi = integrate.dblquad( lambda x,y: dR_det(x,y), 3., 60., lambda x: 3., lambda x: 60. )

print(integ1[0]/denomi[0])
print(integ2[0]/denomi[0])



#plot mass distribution in plane.
d=100
dist_m = np.array([])
sum1 = 0.
sum2 = 0.
for i in range(0,d):
    _m1_ = (i+1) * 1.
    for j in range(0,d):
        _m2_ = (j+1) * 1.
        dist_m = np.append( dist_m, dR_det( _m1_, _m2_ ) )



ax = plt.subplot(111)
im = ax.imshow(dist_m.reshape((d,d)))
# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)

plt.show()






class Metropolis(object):
    
    def __init__(self, func, ndim, proposal_std=1., sample_rate=1):
        self.func = func
        self.ndim = ndim
        self.proposal_std = proposal_std
        self.sample_rate = sample_rate



    def __call__(self, sample_size):
        x = np.full(self.ndim ,30.0)
        samples = []
        for i in range(sample_size * self.sample_rate):
            x_new = np.random.normal(scale=self.proposal_std, size=self.ndim)
            x_new += x
            accept_prob = self.func(x_new) / self.func(x)
            if accept_prob > np.random.uniform():
                x = x_new
            if i % self.sample_rate == 0:
                samples.append(x)
            
            if i % s_s == 0:
                print(sample_size * self.sample_rate - i,'tasks are remaining.')
        assert len(samples) == sample_size
        return np.array(samples)




def func(x):
#    return dR_det( x[0] ,x[1] ,x[2]/100 )
    return dR_det( x[0] ,x[1] )
    


print("R_det sampling...")

sampler = Metropolis(func, 2, proposal_std=15., sample_rate=10)
samples = sampler(sample_size=s_s)

print("mean\n", np.mean(samples, axis=0))
print("covariance\n", np.cov(samples, rowvar=False))

data_gw = np.array([])
data_m1 = np.array([])
data_m2 = np.array([])

print('for example, there are these samples.')


samples_save = np.array([])
"""
for i in range (s_s-10000,s_s):
    samples_save = np.append(samples_save,[samples[i][0],samples[i][1]])
"""
np.save(
    "samples_p=1_alpha=1_lognormal_3over37_10000",  # ファイル名
    samples # 保存したいオブジェクト
)

"""
    box_gw = samples[i]
    print '[',box_gw[0],',',box_gw[1],'],'
    
    data_m1 = np.append(data_m1 ,box_gw[0])
    
    data_m2 = np.append(data_m2 ,box_gw[1])
    
    data_gw = np.append(data_gw , box_gw)
"""
    



plt.scatter(data_m1,data_m2)
plt.xlabel("m1")
plt.ylabel("m2")
plt.grid(True)
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


#for l in range(1, N):
#sum += dR_det( random.uniform( 1.0, 150.0 ) ,random.uniform( 1.0 , 150.0 ) ,random.uniform( 0.01,1.5 ))
            
#print ' denominator is ... '
#print 149 * 149 * 1.49 * sum / N

