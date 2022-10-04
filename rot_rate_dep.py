from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt

NUMBER_OF_STEPS = 1000

MIN_LENGTH = 0.140
MAX_LENGTH = 0.325
LENGTH_STEP = (MAX_LENGTH - MIN_LENGTH) / (NUMBER_OF_STEPS - 1)
na = np.zeros(NUMBER_OF_STEPS)
rel_rot = np.zeros(NUMBER_OF_STEPS)
for i in range(NUMBER_OF_STEPS):
    length = MIN_LENGTH + i * LENGTH_STEP
    tmp_res = rs.Resonator(3, rs.Mirror(0, 180 / np.pi / 1000000, 0.005), rs.Mirror(0.009999, 10, 0.010), rs.Mirror(length, 0, 0.400))
    angle = tmp_res.realign() / 180 * np.pi * 1000000
    if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
        na[i] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
    else:
        na[i] = nan
    rel_rot[i] = angle / na[i]

plt.plot(na, rel_rot)

MIN_RADIUS = 0.301
MAX_RADIUS = 0.400
RADIUS_STEP = (MAX_RADIUS - MIN_RADIUS) / (NUMBER_OF_STEPS - 1)
na = np.zeros(NUMBER_OF_STEPS)
rel_rot = np.zeros(NUMBER_OF_STEPS)
for i in range(NUMBER_OF_STEPS):
    radius = MIN_RADIUS + i * RADIUS_STEP
    tmp_res = rs.Resonator(3, rs.Mirror(0, 180 / np.pi / 1000000, 0.005), rs.Mirror(0.009999, 10, 0.010), rs.Mirror(0.310, 0, radius))
    angle = tmp_res.realign() / 180 * np.pi * 1000000
    if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
        na[i] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
    else:
        na[i] = nan
    rel_rot[i] = angle / na[i]

plt.plot(na, rel_rot)

plt.show()
