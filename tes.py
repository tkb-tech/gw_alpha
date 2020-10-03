print ('hello world!')

import gwdet
import numpy as np
import math
import scipy
from scipy import integrate
from scipy import stats

print('please select sensitivity mode.')
x = input()
if x==1:
    p = gwdet.detectability()
    print('the sensitibity mode is design sensitibity.')

else:
    p = gwdet.detectability(psd = 'aLIGOEarlyHighSensitivityP1200087')
    print('the sensitibity mode is aLIGOEarlyHighSensitivityP1200087.')

#p = gwdet.detectability(psd = 'aLIGOEarlyHighSensitivityP1200087')
#p = gwdet.detectability()

import time

t1 = time.time()

N = 1000

for i in range(0,N):
    print (p(0.05*i,0.05*i,0.1))

t2 = time.time()

print (t2-t1)
