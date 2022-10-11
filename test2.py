from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt

FUND_WAVE = 1064 * 10 ** -7
ANGLE = 0 / 180 * np.pi
R_2 = 10 / 11
R_3 = 200 / 11
chi = 5 / 6
xi = np.min(np.array([2 / np.cos(ANGLE), np.cos(ANGLE) + chi, 1 + chi - R_2 / 2 / np.cos(ANGLE)])) * 0.99999

res1 = rs.Resonator(3, rs.Mirror(0, 0, 1 - R_2 / 2 / chi), rs.Mirror(1, ANGLE, R_2), rs.Mirror(1 + R_3 + R_2 / 2 / (xi - chi), 0, R_3))
res1.waist_scheme(FUND_WAVE)
print(res1.is_g_stable_tangential())
print(res1.is_g_stable_sagittal())
print(res1.waist_search(1064 * 10 ** -9))

res2 = rs.Resonator(2, rs.Mirror(0, 0, 0.005), rs.Mirror(0.00999956, 0, 0.005))
print(res2.waist_search(1064 * 10 ** -9))

