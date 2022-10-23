from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt

FUND_WAVE = 1064 * 10 ** -7
ANGLE = 0 / 180 * np.pi
R_2 = 10 / 11
R_3 = 200 / 11



chi2 = 0.17  #(1 / np.cos(ANGLE) - chi1) * 1.0001    (np.cos(ANGLE) - R_2 / 2) * 0.9999   //  0.17

chi1 = (np.cos(ANGLE) - chi2 * R_2 / (2 * chi2 * R_3 + R_2)) * 0.99999    #5 / 6  //  (1 / np.cos(ANGLE) - chi2) * 1.0001   (np.cos(ANGLE) - chi2 * R_2 / (2 * chi2 * R_3 + R_2)) * 0.99999

res1 = rs.Resonator(3, rs.Mirror(0, 0, 1 - R_2 / 2 / chi1), rs.Mirror(1, ANGLE, R_2), rs.Mirror(1 + R_3 + R_2 / 2 / chi2, 0, R_3))
res1.waist_scheme(FUND_WAVE)
print(res1.waist_search(1064 * 10 ** -7))
plt.show()

