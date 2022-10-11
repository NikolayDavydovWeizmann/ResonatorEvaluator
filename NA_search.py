from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

plt.rcParams.update({'font.size': 30})

NUMBER_OF_LENGTH_STEPS = 100
NUMBER_OF_RADIUS_STEPS = 100

MIN_LENGTH = 0.127
MAX_LENGTH = 0.5
LENGTH_STEP = (MAX_LENGTH - MIN_LENGTH) / (NUMBER_OF_LENGTH_STEPS - 1)

MIN_RADIUS = 0.100
MAX_RADIUS = 0.5
RADIUS_STEP = (MAX_RADIUS - MIN_RADIUS) / (NUMBER_OF_RADIUS_STEPS - 1)

SHORT_ARM_LENGTH = 0.011

na = np.zeros((NUMBER_OF_RADIUS_STEPS, NUMBER_OF_LENGTH_STEPS))
for i in range(NUMBER_OF_LENGTH_STEPS):
    for j in range(NUMBER_OF_RADIUS_STEPS):
        length = MIN_LENGTH + i * LENGTH_STEP
        radius = MIN_RADIUS + j * RADIUS_STEP
        tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 0.005), rs.Mirror(SHORT_ARM_LENGTH, 10 / 180 * np.pi, 0.010), rs.Mirror(length + SHORT_ARM_LENGTH, 0, radius))
        if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
            na[NUMBER_OF_RADIUS_STEPS - 1 - j, i] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
        else:
            na[NUMBER_OF_RADIUS_STEPS - 1 - j, i] = nan

min_na = np.nanmin(na)
max_na = np.nanmax(na)

plt.imshow(na, extent= ([1000 * MIN_LENGTH, 1000 * MAX_LENGTH, 1000 * MIN_RADIUS, 1000 * MAX_RADIUS]), norm= colors.LogNorm(vmin= min_na, vmax= max_na), cmap= cm.get_cmap('viridis'), aspect= 'auto')
plt.xlabel('Length, mm')
plt.ylabel('Radius, mm')
plt.colorbar()
plt.get_current_fig_manager().window.state('zoomed')
plt.show()
