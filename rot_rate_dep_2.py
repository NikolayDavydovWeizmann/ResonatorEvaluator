from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm

plt.rcParams.update({'font.size': 30})
cmap = matplotlib.cm.get_cmap('brg')

NUMBER_OF_LENGTH_STEPS = 100
NUMBER_OF_RADIUS_STEPS = 5

MIN_LENGTH = np.array([0.160, 0.145, 0.130, 0.200, 0.270])
MAX_LENGTH = np.array([0.326, 0.33323, 0.33282, 0.3325, 0.33225])
LENGTH_STEP = (MAX_LENGTH - MIN_LENGTH) / (NUMBER_OF_LENGTH_STEPS - 1)


MIN_RADIUS = 0.340
MAX_RADIUS = 0.600
RADIUS_STEP = (MAX_RADIUS - MIN_RADIUS) / (NUMBER_OF_RADIUS_STEPS - 1)

plt.figure(1)
plt.xlabel('NA')
plt.ylabel(r'$\frac{\Delta \alpha_{mirror}}{\Delta \alpha_{axis} \ / \ NA}$')
plt.grid(visible= True, which= 'both')

na = np.zeros(NUMBER_OF_LENGTH_STEPS)
rel_rot = np.zeros(NUMBER_OF_LENGTH_STEPS)
for i in range(NUMBER_OF_RADIUS_STEPS):
    radius = MIN_RADIUS + i * RADIUS_STEP
    color_of_line = cmap(i / NUMBER_OF_RADIUS_STEPS)
    for j in range(NUMBER_OF_LENGTH_STEPS):
        length = MIN_LENGTH[i] + j * LENGTH_STEP[i]
        tmp_res = rs.Resonator(3, rs.Mirror(0, 10 ** -12, 0.005), rs.Mirror(0.010, 10 / 180 * np.pi, 0.010), rs.Mirror(length, 0, radius))
        angle = tmp_res.realign() * 10 ** 12
        if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
            na[j] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
        else:
            na[j] = nan
        rel_rot[j] = na[j] / angle
    plt.plot(na, rel_rot, label= str(np.around(1000 * radius, 0)) + ' mm', color= color_of_line)
"""
for i in range(NUMBER_OF_LENGTH_STEPS):
    length = MIN_LENGTH + i * LENGTH_STEP
    color_of_line = cmap(i / NUMBER_OF_LENGTH_STEPS)
    for j in range(NUMBER_OF_RADIUS_STEPS):
        radius = MIN_RADIUS[i] + j * RADIUS_STEP[i]
        tmp_res = rs.Resonator(3, rs.Mirror(0, 1 / 1000000, 0.005), rs.Mirror(0.01006, 10 / 180 * np.pi, 0.010), rs.Mirror(length, 0, radius))
        angle = tmp_res.realign() * 1000000
        if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
            na[j] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
        else:
            na[j] = nan
        rel_rot[j] = na[j] / angle
    plt.plot(na, rel_rot, label= str(np.around(1000 * length, 0)) + ' mm', color= color_of_line, linestyle= 'dashed')
"""
plt.legend()
plt.get_current_fig_manager().window.state('zoomed')


plt.figure(2)
plt.xlabel('NA')
plt.ylabel(r'$\frac{\Delta \alpha_{mirror}}{\Delta \alpha_{axis} \ / \ NA}$')
plt.grid(visible= True, which= 'both')

for i in range(NUMBER_OF_RADIUS_STEPS):
    radius = MIN_RADIUS + i * RADIUS_STEP
    color_of_line = cmap(i / NUMBER_OF_RADIUS_STEPS)
    for j in range(NUMBER_OF_LENGTH_STEPS):
        length = MIN_LENGTH[i] + j * LENGTH_STEP[i]
        tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 0.005), rs.Mirror(0.01, 10 / 180 * np.pi, 0.010), rs.Mirror(length, 1 / 1000000, radius))
        angle = tmp_res.realign() * 1000000
        if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
            na[j] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
        else:
            na[j] = nan
        rel_rot[j] = na[j] / angle
    plt.plot(na, rel_rot, label= str(np.around(1000 * radius, 0)) + ' mm', color= color_of_line)
plt.legend()
plt.get_current_fig_manager().window.state('zoomed')


MIN_RADIUS = 0.00450025
MAX_RADIUS = 0.00500365
RADIUS_STEP = (MAX_RADIUS - MIN_RADIUS) / (NUMBER_OF_RADIUS_STEPS - 1)

plt.figure(3)
plt.xlabel('NA')
plt.ylabel(r'$\frac{\Delta \alpha_{mirror}}{\Delta \alpha_{axis} \ / \ NA}$')
plt.grid(visible= True, which= 'both')

for i in range(NUMBER_OF_RADIUS_STEPS):
    color_of_line = cmap(i / NUMBER_OF_RADIUS_STEPS)
    radius = MIN_RADIUS + i * RADIUS_STEP
    MIN_LENGTH[i] = 2 * radius - 0.0000074
    MAX_LENGTH[i] = 2 * radius - 0.0000005
    LENGTH_STEP[i] = (MAX_LENGTH[i] - MIN_LENGTH[i]) / (NUMBER_OF_LENGTH_STEPS - 1)
    for j in range(NUMBER_OF_LENGTH_STEPS):
        length = MIN_LENGTH[i] + j * LENGTH_STEP[i]
        tmp_res = rs.Resonator(2, rs.Mirror(0, 1 / 1000000, radius), rs.Mirror(length, 0, radius))
        angle = tmp_res.realign() * 1000000
        if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
            na[j] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
        else:
            na[j] = nan
        rel_rot[j] = na[j] / angle 
    plt.plot(na, rel_rot, label= str(np.around(1000 * radius, 4)) + ' mm', color= color_of_line)
plt.legend()
plt.get_current_fig_manager().window.state('zoomed')

plt.show()
