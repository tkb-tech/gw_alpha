import gwdet
import numpy as np
import scipy
from scipy import integrate
from scipy import stats
import math
import random
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
import time
import pycbc

"""
# Defaults values
defaults={  'directory' : os.path.dirname(__file__),
            'mcn' : int(1e8),
            'mcbins' : int(1e5),
            'approximant' : 'IMRPhenomD',
            'psd' : 'aLIGOZeroDetHighPower',
            'flow' : 10.,
            'deltaf' : 1./40.,
            'snrthreshold': 8.,
            'massmin' : 1.,
            'massmax' : 100.,
            'zmin' : 1e-4,
            'zmax' : 2.2,
            'mc1d' : int(200)}
"""

#name = 'aLIGOZeroDetHighPower'
name = 'aLIGOLateHighSensitivityP1200087'

delta_f = 1./30.

low_f = 10.

#waveforma_approximant = pycbc.waveform.get_fd_waveform()

psd1 = pycbc.psd.analytical.from_string(name,200000,delta_f,low_f)

x1 = list( np.arange( 0., len(psd1)*delta_f, delta_f ) )

for i in range(len(psd1)):
    psd1[i] =  np.sqrt( psd1[i] )

y1 = psd1

"""
plt.plot(x1,y1)
plt.xlim(5.,10**(4))
plt.yscale('log')
plt.xscale('log')
plt.show()
"""


#name = 'aLIGOEarlyHighSensitivityP1200087'
name = 'aLIGOMidHighSensitivityP1200087'


delta_f = 1./40.

low_f = 10.

#waveforma_approximant = pycbc.waveform.get_fd_waveform()

psd2 = pycbc.psd.analytical.from_string(name,200000,delta_f,low_f)

x2 = list( np.arange( 0., len(psd2)*delta_f, delta_f ) )

for i in range(len(psd2)):
    psd2[i] =  np.sqrt( psd2[i] )

y2 = psd2

plt.plot(x1,y1)
plt.plot(x2,y2)
plt.xlim(5.,10**(4))
plt.yscale('log')
plt.xscale('log')
plt.grid(which='major',color='black',linestyle='-')
plt.grid(which='minor',color='black',linestyle='-')
plt.show()

