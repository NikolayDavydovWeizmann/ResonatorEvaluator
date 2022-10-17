from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt

FUND_WAVE = 1064 * 10 ** -7

ANGLE = 15 / 180 * np.pi

NA = 0.2
z_R = FUND_WAVE / np.pi / NA ** 2

l_1 = 1
L_3 = 100

f = (1 + z_R ** 2) / (1 + ((1 + z_R ** 2) / L_3) + np.sqrt(((1 + z_R ** 2) / L_3)) ** 2 - z_R ** 2)
z_R_1 = z_R * (1 + z_R ** 2) / (((1 + z_R ** 2) / f - 1) ** 2 + z_R ** 2)

r_1 = l_1 + z_R ** 2 / l_1
r_2 = 2 * f
r_3 = L_3 / 2 + z_R_1 ** 2 * 2 / L_3

res1 = rs.Resonator(3, rs.Mirror(0, 0, r_1), rs.Mirror(2 * l_1, 0, r_2), rs.Mirror(2 * l_1 + L_3, 0, r_3))
res1.waist_scheme(FUND_WAVE)
plt.show()
