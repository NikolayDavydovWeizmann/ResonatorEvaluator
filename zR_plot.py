from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt



FUND_WAVE = 1064 * 10 ** -7

NA = 0.2
z_R = FUND_WAVE / np.pi / NA ** 2

f = np.arange(0.01, 1, 0.001)
z_R_1 = z_R * (1 + z_R ** 2) / (((1 + z_R ** 2) / f - 1) ** 2 + z_R ** 2)
l_0_1 = z_R_1 / z_R * ((1 + z_R ** 2) / f - 1)

plt.figure(1)
plt.yscale('log')
plt.plot(f, z_R_1)
plt.figure(2)
plt.yscale('log')
plt.plot(f, l_0_1)

NA = 0.15
z_R = FUND_WAVE / np.pi / NA ** 2

z_R_1 = z_R * (1 + z_R ** 2) / (((1 + z_R ** 2) / f - 1) ** 2 + z_R ** 2)

plt.figure(1)
plt.plot(f, z_R_1)
plt.figure(2)
plt.plot(f, l_0_1)

NA = 0.1
z_R = FUND_WAVE / np.pi / NA ** 2


z_R_1 = z_R * (1 + z_R ** 2) / (((1 + z_R ** 2) / f - 1) ** 2 + z_R ** 2)

plt.figure(1)
plt.plot(f, z_R_1)
plt.figure(2)
plt.plot(f, l_0_1)

plt.show()
