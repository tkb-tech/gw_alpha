print('please input the value of alpha...')

inpt = input()

print('Thanks. The value is',inpt)

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
_alpha_ = 36.0/37.0
s_s = 100
m = [5.0,20.0,30.,50.,100.0]
grid = len(m) - 1
i_a = 1

def poisson(n,v):
    return (v**n) * (np.exp(-v)) / math.factorial(n)

def multi(n,v):
    result = math.factorial(total)
    for i in range(0,grid):
        for j in range(0,grid):
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
samples = [ [ 90.85046241 , 80.47555143 , 67.30139741],
     [ 43.48862358,  49.58930188, 122.13664527],
     [ 30.66299444 , 80.55069217 ,115.66991988],
     [ 13.34255772 , 23.23733218 , 34.76113943],
     [ 40.97871948 , 37.63209779 , 82.22310271],
     [ 21.08710136 , 39.72335203  ,64.31746491],
     [ 26.45203796 ,  9.94536307 , 43.38745374],
     [ 98.43212329 , 87.05852953  ,58.21271271],
     [ 65.53466732 , 29.07740121  ,65.05569387],
     [ 33.51143157 , 12.74901078  ,35.85671442],
     [ 14.07735124 , 62.1998195   ,19.84661972],
     [ 20.1275379  , 14.16375468  , 6.42646762],
     [ 20.19351793 , 28.08857847  ,23.58294904],
     [ 22.90018851 , 30.48298317  ,11.08588433],
     [ 31.92676177 , 57.52545033  ,89.48995244],
     [ 36.81056271 , 14.05255828  ,12.37256808],
     [ 33.33068687 ,  6.30767919  ,10.31475845],
     [ 24.49853201 , 19.72308649  ,30.73562406],
     [ 33.35670567 , 11.50638812  ,56.26369613],
     [ 29.65059889 , 73.10599167  ,52.91328799],
     [ 47.60446631 , 80.71787892  ,85.3541094 ],
     [ 29.08933764 , 49.21055121 ,118.22883792],
     [ 18.74765796 , 39.01671379 , 12.41926334],
     [ 50.9006638  , 40.0414001  , 63.16182618],
     [ 47.18560697 , 48.12459676 , 83.80151042],
     [ 39.69882139 , 37.3594855  , 50.99589711],
     [ 98.11360831 , 29.96077077 , 48.72461022],
     [ 62.00692468 , 51.66796508 ,150.79774747],
     [ 27.23946962 , 36.63544654 , 70.50325421],
     [ 15.66851014 , 32.48700734 , 30.81319371],
     [ 31.95219068 , 22.85867712 , 41.43587276],
     [155.6651991  , 67.1646923  , 65.35171946],
     [ 14.80387223 , 20.16010743 , 26.50806513],
     [ 64.0316951  , 75.91845897 ,108.61166049],
     [ 57.21962774 , 56.1517206  ,102.18269061],
     [ 32.14048797 , 31.81331174 , 63.32441902],
     [ 18.9234606  , 34.89187545 ,  9.80648655],
     [ 62.17172062 , 36.75512998 , 54.98557297],
     [ 19.90354208 , 24.78319873 , 17.17095796],
     [ 55.74504995 , 31.39120965 , 40.35228927],
     [102.93770529 , 63.23564975 , 53.69673856],
     [ 15.44632544 , 14.76428955 , 43.08355419],
     [ 25.02438799 , 31.30282878 , 94.67357993],
     [ 20.99199634 , 13.32028682 , 36.0552651 ],
     [ 56.82851718 , 11.29992839 , 39.56598148],
     [ 37.68373745 , 22.96136394 , 28.93119773],
     [ 53.82492913 , 53.41188783 , 32.54378107],
     [ 37.60345458 , 43.82699461 , 26.85852819],
     [ 29.83641447 ,172.46691215 , 44.57932349],
     [ 55.24164608 , 12.0102583  , 42.99938948],
     [ 40.16889857 , 30.45605693 , 13.45772885],
     [ 19.67581485 , 52.88706489 , 41.39821449],
     [ 47.66123754 , 27.66399727 , 29.4327272 ],
     [ 94.12334656 , 21.00368769 , 30.49390945],
     [ 32.47338894 , 17.54902593 , 22.71276402],
     [ 18.2181427  , 12.10588177 , 15.41243416],
     [ 66.00502248 , 22.26660986 , 33.02605891],
     [115.67139157 , 34.69416316 , 38.35080674],
     [ 51.96309466 , 73.67360861 , 97.18636348],
     [ 41.73994363 , 37.44302329 , 78.102715  ],
     [ 79.61112298 , 27.46760912 , 61.2332272 ],
     [ 35.17578514 , 41.62272737 , 15.77416587],
     [ 35.66097929 , 53.26138177 , 19.15578673],
     [ 29.30516489 , 87.65306593 , 51.80294363],
     [ 81.04505131 , 26.77682661 , 30.87855121],
     [ 44.82942546 , 14.32319132 , 22.36743338],
     [ 81.38959539 , 15.87239752 , 39.37027637],
     [ 68.85590803 , 56.59405562 , 58.66269847],
     [ 27.00436357 , 58.48644861 , 14.46566448],
     [ 49.09229999 , 62.38521341 , 51.74848044],
     [ 46.46615553 , 48.42483329 , 80.29130784],
     [ 16.12253796 , 24.15726615 , 43.90758116],
     [ 44.10992952 , 35.87036489 , 93.78241073],
     [ 37.65423597 ,147.29344157 , 84.35069213],
     [ 21.96352242 , 92.73827003 , 81.54278129],
     [ 83.2010068  , 35.33678601 , 67.17409346],
     [ 45.41647001 , 62.72409169 , 89.3586355 ],
     [ 25.7093488  , 63.56986196 , 37.55573623],
     [ 46.63277362 , 36.94107285 , 21.35154077],
     [ 22.11430055 , 48.37710664 , 23.93742739],
     [ 27.94271308 , 27.13623448 , 23.48794111],
     [ 50.31286968 , 40.84369493 ,102.36884008],
     [ 61.46612107 , 45.40351282 , 74.02113032],
     [ 63.91848402 , 25.62120877 , 45.00765724],
     [ 16.76476911 , 37.29591916 , 80.25314642],
     [ 18.15242032 , 53.7077267  , 25.89320319],
     [ 26.23169925 , 54.05228471 , 27.71948339],
     [ 35.00731374 , 29.79370649 , 16.49110217],
     [ 33.46812145 , 48.93298116 , 85.03579367],
     [ 38.59441877 , 52.14207688 , 15.89093589],
     [ 69.44939668 ,100.31831765 ,108.45980414],
     [ 43.75106785 , 40.94086936 , 76.78891143],
     [ 64.26942635 , 33.35382064 , 81.1237083 ],
     [ 68.35198568 , 27.52883142 , 91.64840082],
     [ 45.76820763 , 33.34218004 , 44.81113355],
     [ 76.59388317 , 72.48843118 , 97.48432955],
     [ 60.66792861 , 26.82526267 , 45.05999242],
     [ 34.94864317 , 22.57308379 , 36.54895968],
     [ 61.12308118 , 21.91134655 , 19.19825399],
     [ 22.20698363 , 78.62613243 , 53.16239133] ]

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

print samples
mean_redshift = np.mean(samples, axis=0)
print "mean\n", mean_redshift
print "covariance\n", np.cov(samples, rowvar=False)



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
Const_in_R_th = 1.0
"""
def f(m):
    return ( 1/ ( np.sqrt(2.0*np.pi) * sigma * m ) ) * math.exp( -( math.log(m/m_c) )**2.0 / ( 2.0 * sigma**2.0 ) )
"""

def f(m):
    if 5.0 < m and m < 60.0:
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

#    p = gwdet.detectability()
    def p(m1,m2,z):
        return 1.0

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
    z_i = 0.00
    z_f = 2.0
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

    p_n_alpha = 1.0
    for i in range(0,grid):
        for j in range(0,grid):
            p_n_alpha = p_n_alpha * poisson(n[i][j],v[i][j])
    print 'the p_n_alpha by poisson is'
    print p_n_alpha

    baki = multi(n,v)
    print 'the p_n_alpha by multinomial is'
    print baki
    probability.append([_alpha_,p_n_alpha])
    probability_multinomial.append([_alpha_,baki])
    print ' '

t2 = time.time()

print res_kai
print probability
print probability_multinomial
print 'the mass function is ;flat'
print 'masses that splits the plane are'
print m


print 'the total time is'
print (t2-t1)/60.0
print 'minutes'
