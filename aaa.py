import matplotlib.pyplot as plt
import numpy as np


class Metropolis(object):
    
    def __init__(self, func, ndim, proposal_std=1., sample_rate=1):
        self.func = func
        self.ndim = ndim
        self.proposal_std = proposal_std
        self.sample_rate = sample_rate
    
    def __call__(self, sample_size):
        x = np.zeros(self.ndim)
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
    return np.exp(-0.5 * np.sum(x ** 2, axis=-1) / 5.)
    
print "\ntwo dimensional"
    
sampler = Metropolis(func, 2, proposal_std=2., sample_rate=10)
samples = sampler(sample_size=100)
    
print "mean\n", np.mean(samples, axis=0)
print "covariance\n", np.cov(samples, rowvar=False)
    
x, y = np.meshgrid(
    np.linspace(-10, 10, 100), np.linspace(-10, 10, 100))
z = func(np.array([x, y]).reshape(2, -1).T).reshape(100, 100)
plt.contour(x, y, z)
plt.scatter(samples[:, 0], samples[:, 1])
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.show()

