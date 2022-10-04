from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt

NUM_OF_LENGTH_STEPS = 100
MIN_LENGTH = 0.200
MAX_LENGTH = 0.350
LENGTH_STEP = (MAX_LENGTH - MIN_LENGTH) / (NUM_OF_LENGTH_STEPS - 1)

NUM_OF_RADIUS_STEPS = 100
MIN_RADIUS = 0.100
MAX_RADIUS = 0.700
RADIUS_STEP = (MAX_RADIUS - MIN_RADIUS) / (NUM_OF_RADIUS_STEPS - 1)

na1 = np.zeros(NUM_OF_RADIUS_STEPS * NUM_OF_LENGTH_STEPS)
na2 = np.zeros(NUM_OF_RADIUS_STEPS * NUM_OF_LENGTH_STEPS)
length_of_II_arm = np.zeros(NUM_OF_RADIUS_STEPS * NUM_OF_LENGTH_STEPS)
radius_of_term_mirror = np.zeros(NUM_OF_RADIUS_STEPS * NUM_OF_LENGTH_STEPS)

for i in range(NUM_OF_LENGTH_STEPS):
    for j in range(NUM_OF_RADIUS_STEPS):
        indx = i * NUM_OF_RADIUS_STEPS + j
        offset = i * LENGTH_STEP
        length_of_II_arm[indx] = MIN_LENGTH + offset
        radius_of_term_mirror[indx] = MIN_RADIUS + j * RADIUS_STEP
        tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 0.005), rs.Mirror(0.009999, 10, 0.010), rs.Mirror(MIN_LENGTH + offset, 0, radius_of_term_mirror[indx]))
        if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
            na1[indx] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
            na2[indx] = tmp_res.waist_search(1064 * 10 ** -9)[0, 1, 2]
        else:
            na1[indx] = nan
            na2[indx] = nan

plt.figure(1)
ax = plt.axes(projection= '3d')
ax.set_title('Numerical aperture (NA) (I arm)')
ax.set_xlabel('Length of II arm, mm')
ax.set_ylabel('Radius of terminal mirror, mm')
ax.set_zlabel('NA')
ax.scatter3D(1000 * length_of_II_arm, 1000 * radius_of_term_mirror, na1, s=1)
#ax.scatter3D(1000 * length_of_II_arm, 1000 * radius_of_term_mirror, na2, s=1)

plt.show()
