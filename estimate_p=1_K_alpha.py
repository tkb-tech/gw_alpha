import gwdet
import numpy as np
import scipy
from scipy import integrate
import math
import random
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time

H_0 = 70.0 #[(km/s)/Mpc]
c = 299792.0 #[km/s]
o_M = 0.3
o_k = 0.0
o_l = 0.7
s_s = 10000

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


#p = gwdet.detectability()

error = []

result = []

for theta in range(0,20):
    alpha = 1. + 0.02*theta - 0.2

    def dR_det(m_1,m_2):
        if m_1 <= 0.0:
            return 0.0
        elif m_2 <= 0.0:
            return 0.0
        else:
            return ( m_1 + m_2 )**( alpha ) * h( m_1 ) * h( m_2 )

    def dR_dmt_dq(m_t,q):
        return m_t/(1+q)**2*dR_det(m_t*q/(1+q),m_t/(1+q))




    #calculate the average of m1,m2 in dR_det.

    integ1 = integrate.dblquad( lambda x,y: x*dR_det(x,y) , 0.1 , 150. , lambda x: 0.1 , lambda x:150. )

    integ2 = integrate.dblquad( lambda x,y: (x+y)*dR_det(x,y) , 0.1 , 150. , lambda x: 0.1 , lambda x:150. )

    denomi = integrate.dblquad( lambda x,y: dR_det(x,y), 0.1 , 150., lambda x: 0.1 , lambda x: 150. )

    print(integ1[0]/denomi[0])
    print(integ2[0]/denomi[0])




    #plot mass distribution in plane.
    n_m=700
    n_q=700
#    n_m = 500
#    n_q = 500
#    n_m = 2000
#    n_q = 1500
    
    m_t_min = 10.
    m_t_max = 100.
    dm_t = (m_t_max - m_t_min)/n_m

    q_min = 0.1
    q_max = 0.95
    dq = (q_max - q_min)/n_q
    sum_int = 0.
    dist_m = np.array([])
    mav = 0.
    sum_logmt = 0.
    logmt = np.array([])
    t1 = time.time()
    
    for i in range(0,n_m):
        m_t = dm_t * (i) + m_t_min
        for j in range(0,n_q):
            q = dq * (j) + q_min
            valu = (dR_dmt_dq( m_t , q ) + dR_dmt_dq( m_t + dm_t , q ) + dR_dmt_dq( m_t , q + dq ) + dR_dmt_dq( m_t + dm_t , q + dq ) + dR_dmt_dq( m_t + dm_t/2. , q + dq/2. ) )/ 5.
            logu = valu * math.log(m_t)
            logmt = np.append( logmt, logu )
            sum_logmt = sum_logmt + logu
            dist_m = np.append( dist_m, valu )
            sum_int = sum_int + valu
            
    sum_logmt = sum_logmt/sum_int

    for i in range(0,dist_m.shape[0]):
        dist_m[i] = dist_m[i]/sum_int
        logmt[i] = logmt[i]/sum_int

    t2 = time.time()

    print('the time to calculate is')
    print(t2-t1)

    logmatrix = logmt.reshape((n_m,n_q))
    matrix = dist_m.reshape((n_m,n_q))

    np.save(
        "p=1_lognormal_3_37_10.<mt<100._0.1<q<0.95_alpha="+str(alpha)+"_"+str(n_m)+"_"+str(n_q),  # ファイル名
        matrix # 保存したいオブジェクト
    )
    np.save(
        "p=1_lognormal_3_37_10.<mt<100._0.1<q<0.95_alpha="+str(alpha)+"_"+str(n_m)+"_"+str(n_q)+"_log",  # ファイル名
        logmatrix # 保存したいオブジェクト
    )
    

    """
    t1 = time.time()
    matrix = np.load('10.<mt<100._0.1<q<0.95.npy')
    logmatrix = np.load('10.<mt<100._0.1<q<0.95_log.npy')
    t2 = time.time()

    sum_logmt = np.sum(logmatrix)

    print('the loading time is')
    print(t2-t1)
    """

    """
    ax = plt.subplot(111)

    im = ax.imshow(matrix)
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    plt.colorbar(im, cax=cax)

    plt.show()
    """



    presamples =np.load('samples_p=1_alpha=1_lognormal_3over37_10000.npy')
    print(presamples.shape)
    samples = np.zeros((s_s,2))
    for i in range(len(presamples)-s_s,len(presamples)):
        samples[i-len(presamples)+s_s][0] = presamples[i][0]
        samples[i-len(presamples)+s_s][1] = presamples[i][1]
    print(samples.shape)




    delta_chi = []
    plot_1 = []
    plot_2 = []

    for tkb in range(2,8):
            print(' ')
            nm=tkb
            nq=tkb

            split_m = []

            split_q = []

            i_ = 0
            sumq = []
            while len(split_m) < nm - 1:
                x = 0.
                for j in range(0,n_q): x = x + matrix[i_][j]
                sumq.append(x)
                
                if sum(sumq)>1./nm:
                    if 1./nm - (sum(sumq) - x) >= sum(sumq) - 1./nm:
                        split_m.append(i_)
                    else:
                        split_m.append(i_-1)
                    sumq = []
                i_ = i_ + 1


            split_m.append(n_m-1)
            split_m.insert(0,0)

    #        print(split_m)



            for k in range(0,nm):
                split_q.append([])
                j_ = 0
                summ = []
                while len(split_q[k]) < nq - 1:
                    x = 0.
                    for i in range(split_m[k],split_m[k+1]): x = x + matrix[i][j_]
                    summ.append(x)
                    if sum(summ)>1./(nm*nq):
                        if 1./(nm*nq) - (sum(summ) - x) >= sum(summ) - 1./(nm*nq):
                            split_q[k].append(j_)
                        else:
                            split_q[k].append(j_-1)
                        summ = []
                    j_ = j_ + 1

                split_q[k].append(n_q-1)
                split_q[k].insert(0,0)
    #        print(split_q)


            x = 0.
            w = 0.
            grid_expectation_value = np.zeros((nm,nq))
            sum_log_i = []
            for k in range(0,nm):
                for l in range(0,nq):
                    x = 0.
#                    y = 0.
                    
                    for i in range(split_m[k],split_m[k+1]):
                        for j in range(split_q[k][l],split_q[k][l+1]):
                            x = x + matrix[i][j]
#                            y = y + logmatrix[i][j]
                    if x * nm * nq < 0.8:
                        print('error is...')
                        print(x,x-1./(nm*nq),k,l)
                        error.append([nm*nq,x,x-1./(nm*nq),k,l])
                        
#                    sum_log_i.append(y/x)
                    grid_expectation_value[k][l] = x
                    w = w + x
            print('sum of grid expwctation value...')
            print(w)
            
            
            for i in range(0,nm):
                for j in range(0,nq):
                    grid_expectation_value[i][j] = s_s * grid_expectation_value[i][j] / w
            
            
            
            
            """
            print('log expectation value for each grid is')
            print(sum_log_i)
            print('sum of it is')
            print(sum(sum_log_i))

            print('expectation value of log(mt) is')
            print(sum_logmt)

            print('numeriator of d\chi is')
            print(sum(sum_log_i)-nm*nq*sum_logmt)
            print((sum(sum_log_i)-nm*nq*sum_logmt)/np.sqrt(2*(nm*nq-1)))
            print(' ')
            """
            delta_chi.append([nm,nq,(sum(sum_log_i)-nm*nq*sum_logmt)/np.sqrt(2*(nm*nq-1))])
            
            
            
            
            plot_1.append(nm)
            plot_2.append((sum(sum_log_i)-nm*nq*sum_logmt)/np.sqrt(2*(nm*nq-1)))
            
            
                
            n = np.zeros((nm,nq))
            total = 0.0
            sum_v = 0.0

            for i in range(0,nm):
                for j in range(0,nq):
                    for k in range(0,s_s):
                        m_1 = samples[k][0]
                        m_2 = samples[k][1]
#                        z = samples[k][2]
                        m_total = m_1 + m_2
                        if m_1 < m_2:
                            q = m_1/m_2
                        else:
                            q = m_2/m_1
                        
                        if dm_t * split_m[i] + m_t_min < m_total and m_total < dm_t * split_m[i+1] + m_t_min and dq * split_q[i][j] + q_min < q and q < dq * split_q[i][j+1] + q_min:
                            n[i][j] = n[i][j] + 1.0
                            total = total + 1.0
            
            
            for i in range(0,nm):
                for j in range(0,nq):
                    grid_expectation_value[i][j] = grid_expectation_value[i][j] * (total/s_s)
                    
            print('total is')
            print(total)
            print('sum of grid_expectation_value is')
            print(np.sum(grid_expectation_value))
            
            kai_n = 0.0
            for i in range(0,nm):
                for j in range(i,nq):
                    kai_n = kai_n + (n[i][j]-grid_expectation_value[i][j])**2.0/grid_expectation_value[i][j]
            print('when the K = ',nm*nq)
            print('the chi value is')
            print(kai_n)
            print(' ')
            result.append([nm*nq,alpha,kai_n])


print(result)

hoge = np.array(result)

np.save('11_3_result',hoge)

"""
    #nm = 3
    #nq = 3
    #nm = 4
    #nq = 4
    #nm = 5
    #nq = 5
    #nm = 6
    #nq = 6
    nm=7
    nq=7

    split_m = []

    split_q = []

    i_ = 0
    sumq = []
    while len(split_m) < nm - 1:
        x = 0.
        for j in range(0,n_q): x = x + matrix[i_][j]
        sumq.append(x)
        
        if sum(sumq)>1./nm:
            if 1./nm - (sum(sumq) - x) >= sum(sumq) - 1./nm:
                split_m.append(i_)
            else:
                split_m.append(i_-1)
            sumq = []
        i_ = i_ + 1


    split_m.append(n_m-1)
    split_m.insert(0,0)

    print(split_m)



    for k in range(0,nm):
        split_q.append([])
        j_ = 0
        summ = []
        while len(split_q[k]) < nq - 1:
            x = 0.
            for i in range(split_m[k],split_m[k+1]): x = x + matrix[i][j_]
            summ.append(x)
            if sum(summ)>1./(nm*nq):
                if 1./(nm*nq) - (sum(summ) - x) >= sum(summ) - 1./(nm*nq):
                    split_q[k].append(j_)
                else:
                    split_q[k].append(j_-1)
                summ = []
            j_ = j_ + 1

        split_q[k].append(n_q-1)
        split_q[k].insert(0,0)
    print(split_q)


    x = 0.

    sum_log_i = []
    for k in range(0,nm):
        for l in range(0,nq):
            x = 0.
            y = 0.
            
            for i in range(split_m[k],split_m[k+1]):
                for j in range(split_q[k][l],split_q[k][l+1]):
                    x = x + matrix[i][j]
                    y = y + logmatrix[i][j]
            sum_log_i.append(y/x)
            print(x)
    print('log expectation value for each grid is')
    print(sum_log_i)
    print('sum of it is')
    print(sum(sum_log_i))

    print('expectation value of m_t is')
    print(mav/sum_int)
    print('expectation value of log(mt) is')
    print(sum_logmt/sum_int)


    uhei = 0.

    for i in range(0,n_m):
        for j in range(0,n_q):
            uhei = uhei + logmatrix[i][j]

    print('numeriator of d\chi is')
    print(sum(sum_log_i)-nm*nq*sum_logmt/sum_int)
    print(sum(sum_log_i)-nm*nq*uhei)
"""




#次にやるべきことは、グリッドを決める作業ね。p=1の場合だけど...
#p≠1の場合だと、そもそもmatrixを求めるのに莫大な計算時間を要する可能性もあるね。。どうしようかね。まぁそこは仕方ないってことで。
#後やるべきなのは、p=1の時の分布関数とp≠1の時の分布関数にどの程度の差があるのかは知っておきたいよね。もし差がないなら最適なKの評価はp=1の場合でやれば良さげだしね。（もちろん質量関数によるか。。。）
