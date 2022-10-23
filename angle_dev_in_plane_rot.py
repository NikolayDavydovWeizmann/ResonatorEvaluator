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


rot = np.empty(L_3.size)
for i in range(L_3.size):
    res1 = rs.Resonator(3, rs.Mirror(0, 0, r_1, in_plane_angle_deviation=delta),
                        rs.Mirror(2 * l_1, 0, r_2[i]),
                        rs.Mirror(2*l_1 + L_3[i], 0, r_3[i]))
    rot[i] = res1.realign()

res2 = rs.Resonator(2, rs.Mirror(0, 0, r_1, in_plane_angle_deviation=delta),
                    rs.Mirror(2 * l_1, 0, r_1))
rot2 = res2.realign()

plt.figure(1)
plt.plot(L_3 / 2, rot2 / rot)

plt.figure(2)
plt.plot(L_3 / 2, delta / rot)
print(delta / rot2)
# plt.plot([L_3[0] / 2, L_3[-1] / 2], [delta / rot2, delta / rot2])
plt.show()
