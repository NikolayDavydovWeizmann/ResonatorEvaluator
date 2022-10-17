from math import nan
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

plt.rcParams.update({'font.size': 30})

ANGLE = 0 / 180 * np.pi
FUND_WAVE = 1064 * 10 ** -7
F = 1 / 2 * np.cos(ANGLE)

NUMBER_OF_L1_STEPS = 10000
NUMBER_OF_L2_STEPS = 10000

MIN_L1 = 1.0001 * F
MAX_L1 = 0.6
L1_STEP = (MAX_L1 - MIN_L1) / (NUMBER_OF_L1_STEPS - 1)

MIN_L2 = 1.0001 * F
MAX_L2 = 5
L2_STEP = (MAX_L2 - MIN_L2) / (NUMBER_OF_L2_STEPS - 1)

na = np.zeros([NUMBER_OF_L2_STEPS, NUMBER_OF_L1_STEPS])

for i in range(NUMBER_OF_L1_STEPS):
    length1 = MIN_L1 + i * L1_STEP
    for j in range(NUMBER_OF_L2_STEPS):
        length2 = MIN_L2 + j * L2_STEP
        tmp_waist2 = FUND_WAVE / np.pi * np.sqrt((length1 + length2 - length1 * length2 / F) * F * (length1 - F) / (length2 - F))
        na[NUMBER_OF_L2_STEPS - 1 - j, i] = FUND_WAVE / np.pi / np.sqrt(tmp_waist2)

min_na = np.nanmin(na)
max_na = np.nanmax(na)

plt.title(r"$NA_{tg}$")
plt.imshow(na, extent= ([MIN_L1, MAX_L1, MIN_L2, MAX_L2]), norm= colors.LogNorm(vmin= min_na, vmax= max_na), cmap= cm.get_cmap('viridis'), aspect= 'auto')
plt.xlabel(r'$l_1$')
plt.ylabel(r'$l_2$')
plt.colorbar()
plt.get_current_fig_manager().window.state('zoomed')

plt.show()