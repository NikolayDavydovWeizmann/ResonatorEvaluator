from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt


FUND_WAVE = 1064 * 10 ** -7
NA = 0.045
l_2 = 25

r_1  = 1 / 2 * (1  + 4 * FUND_WAVE ** 2 / (np.pi ** 2 * NA ** 4))


NA_2 = NA / l_2

r_3 = l_2 / 2 * (1  + 4 * FUND_WAVE ** 2 / (np.pi ** 2 * NA_2 ** 4 * l_2 ** 2))


r_2 = 2 * r_1 * r_3 / (r_1 + r_3)


res1 = rs.Resonator(2, rs.Mirror(0, 0, r_1), rs.Mirror(1, 0, r_1))
print(res1.waist_search(FUND_WAVE))
res1.waist_scheme(FUND_WAVE)

res2 = rs.Resonator(3, rs.Mirror(0, 0, r_1), rs.Mirror(1, 2 / 180 * np.pi, r_2), rs.Mirror(1 + l_2, 0, r_3))
print(res2.waist_search(FUND_WAVE))
res2.waist_scheme(FUND_WAVE)

plt.show()

res1 = rs.Resonator(2, rs.Mirror(0, 10 ** -12, r_1), rs.Mirror(1, 0, r_1))


res2 = rs.Resonator(3, rs.Mirror(0, 10 ** -12, r_1), rs.Mirror(1,  2 / 180 * np.pi, r_2), rs.Mirror(1 + l_2, 0, r_3))
print(res1.realign() / res2.realign())
