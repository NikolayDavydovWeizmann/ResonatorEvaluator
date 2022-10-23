import resonator as rs
import numpy as np
import matplotlib.pyplot as plt


FUND_WAVE = 1064 * 10 ** -7
delta = 10 ** -10

"""
ANGLE = 15 / 180 * np.pi
"""

NA = 0.2
z_R = FUND_WAVE / np.pi / NA**2

l_1 = 1
L_3 = np.arange(2, 1000, 2)

f = (1 + z_R**2) / (1 + (1 + z_R**2)/L_3
                    + np.sqrt(((1 + z_R**2) / L_3)**2 - z_R**2))
z_R_1 = z_R * (1 + z_R**2) / (((1 + z_R**2)/f - 1)**2 + z_R**2)

r_1 = l_1 + z_R**2 / l_1
r_2 = 2 * f
r_3 = L_3 / 2 + z_R_1**2 * 2 / L_3


na = np.empty(L_3.size)
for i in range(L_3.size):
    res1 = rs.Resonator(3, rs.Mirror(0, 0, r_1, in_plane_angle_deviation=delta),
                        rs.Mirror(2 * l_1, 0, r_2[i]),
                        rs.Mirror(2*l_1 + L_3[i], 0, r_3[i]))
    na_tmp = res1.waist_search(FUND_WAVE)[0, 0, 2]
    res1.realign()
    na[i] = res1.waist_search(FUND_WAVE)[0, 0, 2]/na_tmp - 1

res2 = rs.Resonator(2, rs.Mirror(0, 0, r_1, in_plane_angle_deviation=delta),
                    rs.Mirror(2 * l_1, 0, r_1))
na2_init = res2.waist_search(FUND_WAVE)[0, 0, 2]
res2.realign()
na2 = res2.waist_search(FUND_WAVE)[0, 0, 2]/na2_init - 1
print(delta / na2)

plt.figure(1)
plt.yscale('log')
plt.plot(L_3 / 2 / l_1, np.abs(na2 / na))

plt.figure(2)
plt.yscale('log')
plt.plot(L_3 / 2 / l_1, np.abs(delta / na))
plt.plot([L_3[0] / 2, L_3[-1] / 2], [delta / np.abs(na2), delta / np.abs(na2)])
plt.show()
