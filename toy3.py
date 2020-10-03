import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

mu=10
sigma=1
data=np.random.normal(mu,sigma,100)

x = np.linspace(0,20,num = 1000)

kde = stats.gaussian_kde(data)
density = kde(x)

print kde(10)
print 'hoge'
print data
print density


plt.plot(x,kde(x))
plt.show()

#fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
#x, y, z = data
#ax.scatter(x, y, z, c=density)
#plt.show()



