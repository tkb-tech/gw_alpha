import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

data = np.random.normal(10,5,300)

print(data)
#[7.4,
# 7.2,
# 8.0,
# 8.5,
# ...
# 6.2,
# 3.8]

plt.hist(data, bins=31, range=(0, 31))
plt.show()


kde = gaussian_kde(data)


estimates = kde(np.linspace(0, 31, num=32))
print(estimates)
# [0.00409969 0.00757112 0.0128535  0.01999006 0.02838034 0.03670379
#  0.04325381 0.04661941 0.04635048 0.04318489 0.03868037 0.03447436
#  0.03161084 0.03028665 0.03006976 0.03036525 0.03082537 0.03149862
#  0.03267868 0.03457151 0.03699876 0.03934671 0.04082662 0.04088853
#  0.03948617 0.03698289 0.03378769 0.03005188 0.02569315 0.02069462
#  0.01536906 0.01032863]


plt.plot(kde(np.linspace(0, 31, num=32)))
plt.xlabel('temerature [C]')
plt.ylabel('kde')
plt.show()
