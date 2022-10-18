import numpy as np
import matplotlib.pyplot as plt

FUND_WAVE = 1064 * 10 ** -7

NA = 0.05
z_R = FUND_WAVE / np.pi / NA ** 2

l_1 = 1
L_3 = np.arange(2, 100, 2)


f = (1 + z_R ** 2) / (1 + ((1 + z_R ** 2) / L_3) + np.sqrt(((1 + z_R ** 2) / L_3)) ** 2 - z_R ** 2)

z_R_1 = z_R * (1 + z_R ** 2) / (((1 + z_R ** 2) / f - 1) ** 2 + z_R ** 2)
l_0_1 = z_R_1 / z_R * ((1 + z_R ** 2) / f - 1)


xi_1 = 1 - z_R ** 2 / l_1
xi_2 = l_0_1 - 2 * z_R_1 ** 2 / L_3
r_1 = l_1 + z_R ** 2 / l_1
r_3 = L_3 / 2 + z_R_1 ** 2 * 2 / L_3

acc = (f * xi_2 - xi_1 * (xi_2 - f)) / (r_1 * (xi_2 - f))
acc3 = (f * xi_2 - xi_1 * (xi_2 - f)) / r_3 / f


acc_2mirror = 2 * (1 - 1 / r_1)

plt.figure(1)
plt.plot(L_3 / 2, acc / acc_2mirror)
plt.figure(2)
plt.plot(L_3 / 2, acc3 / acc_2mirror)
plt.show()
