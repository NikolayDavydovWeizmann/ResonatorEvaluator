from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})

NUMBER_OF_LENGTH_STEPS = 10
NUMBER_OF_RADIUS_STEPS = 100

MIN_LENGTH = 0.170
MAX_LENGTH = 0.325
LENGTH_STEP = (MAX_LENGTH - MIN_LENGTH) / (NUMBER_OF_LENGTH_STEPS - 1)

MIN_RADIUS = np.array([0.320, 0.289, 0.279, 0.272, 0.270, 0.274, 0.280, 0.288, 0.298, 0.310077])
MAX_RADIUS = np.array([0.530, 0.530, 0.530, 0.530, 0.560, 0.570, 0.580, 0.570, 0.580, 0.400])
RADIUS_STEP = (MAX_RADIUS - MIN_RADIUS) / (NUMBER_OF_RADIUS_STEPS - 1)

plt.figure(1)
plt.xlabel('NA')
plt.ylabel(r'$\frac{\Delta \alpha_{axis}}{NA \cdot \Delta x}, \mu m^{-1}$')

na = np.zeros(NUMBER_OF_RADIUS_STEPS)
rel_rot = np.zeros(NUMBER_OF_RADIUS_STEPS)
for i in range(NUMBER_OF_LENGTH_STEPS):
    length = MIN_LENGTH + i * LENGTH_STEP
    for j in range(NUMBER_OF_RADIUS_STEPS):
        radius = MIN_RADIUS[i] + j * RADIUS_STEP[i]
        tmp_res = rs.Resonator(3, rs.Mirror(0, np.arcsin(10 ** -15 / 0.005), 0.005), rs.Mirror(0.009999, 10 / 180 * np.pi, 0.010), rs.Mirror(length, 0, radius))
        angle = tmp_res.realign() / (10 ** -15)
        if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
            na[j] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
        else:
            na[j] = nan
        rel_rot[j] = angle / na[j]
    plt.plot(na, rel_rot / 1000000, label= str(np.around(1000 * length, 0)) + ' mm')
plt.legend()


plt.figure(2)
plt.xlabel('NA')
plt.ylabel(r'$\frac{\Delta \alpha_{axis}}{NA \cdot \Delta x}, \mu m^{-1}$')

for i in range(NUMBER_OF_LENGTH_STEPS):
    length = MIN_LENGTH + i * LENGTH_STEP
    for j in range(NUMBER_OF_RADIUS_STEPS):
        radius = MIN_RADIUS[i] + j * RADIUS_STEP[i]
        tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 0.005), rs.Mirror(0.009999, 10 / 180 * np.pi, 0.010), rs.Mirror(length, np.arcsin(10 ** -12 / radius), radius))
        angle = tmp_res.realign() / (10 ** -12)
        if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
            na[j] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
        else:
            na[j] = nan
        rel_rot[j] = angle / na[j]
    plt.plot(na, rel_rot / 1000000, label= str(np.around(1000 * length, 0)) + ' mm')
plt.legend()


MIN_LENGTH = 0.009
MAX_LENGTH = 0.0099999
LENGTH_STEP = (MAX_LENGTH - MIN_LENGTH) / (NUMBER_OF_LENGTH_STEPS - 1)

plt.figure(3)
plt.xlabel('NA')
plt.ylabel(r'$\frac{\Delta \alpha_{axis}}{NA \cdot \Delta x}, \mu m^{-1}$')

for i in range(NUMBER_OF_LENGTH_STEPS):
    length = MIN_LENGTH + i * LENGTH_STEP
    MIN_RADIUS[i] = (length + 0.0000074) / 2
    MAX_RADIUS[i] = (length + 0.0000005) / 2
    RADIUS_STEP[i] = (MAX_RADIUS[i] - MIN_RADIUS[i]) / (NUMBER_OF_RADIUS_STEPS - 1)
    for j in range(NUMBER_OF_RADIUS_STEPS):
        radius = MIN_RADIUS[i] + j * RADIUS_STEP[i]
        tmp_res = rs.Resonator(2, rs.Mirror(0, np.arcsin(10 ** -15 / radius), radius), rs.Mirror(length, 0, radius))
        angle = tmp_res.realign() / 10 ** -15
        if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
            na[j] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
        else:
            na[j] = nan
        rel_rot[j] = angle / na[j]
    plt.plot(na, rel_rot / 1000000, label= str(np.around(1000 * length, 4)) + ' mm')
plt.legend()

plt.show()