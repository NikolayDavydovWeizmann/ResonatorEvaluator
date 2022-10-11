from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt


NUM_OF_LENGTH_STEPS = 100
NUM_OF_RADIUS_STEPS = 100
LENGTH_STEP = 0.0005
RADIUS_STEP = 0.003

ld = np.zeros(NUM_OF_LENGTH_STEPS * NUM_OF_RADIUS_STEPS)
radius = np.zeros(NUM_OF_LENGTH_STEPS * NUM_OF_RADIUS_STEPS)
na1 = np.zeros(NUM_OF_LENGTH_STEPS * NUM_OF_RADIUS_STEPS)
na2 = np.zeros(NUM_OF_LENGTH_STEPS * NUM_OF_RADIUS_STEPS)
for i in range(NUM_OF_LENGTH_STEPS):
    for j in range(NUM_OF_RADIUS_STEPS):
        indx = i * NUM_OF_RADIUS_STEPS + j
        ld[indx] = (i - np.around(NUM_OF_LENGTH_STEPS / 2)) * LENGTH_STEP
        radius[indx] = (j - np.around(NUM_OF_RADIUS_STEPS / 2)) * RADIUS_STEP
        tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 0.005), rs.Mirror(0.0099999, 10 / 180 * np.pi, 0.010), rs.Mirror(0.300 + ld[indx], 0, 0.595 + radius[indx]))
        if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
            na1[indx] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
            na2[indx] = tmp_res.waist_search(1064 * 10 ** -9)[1, 0, 2]
        else:
            na1[indx] = nan
            na2[indx] = nan
        
plt.figure(1)
ax = plt.axes(projection= '3d')
ax.set_title('Numerical aperture (NA) (I shoulder)')
ax.set_xlabel('Length deviation, mm')
ax.set_ylabel('Radius deviation, mm')
ax.set_zlabel('NA')
ax.scatter3D(1000 * ld, 1000 * radius, na1, s=1)

plt.figure(2)
ax = plt.axes(projection= '3d')
ax.set_title('Numerical aperture (NA) (II shoulder)')
ax.set_xlabel('Length deviation, mm')
ax.set_ylabel('Radius deviation, mm')
ax.set_zlabel('NA')
ax.scatter3D(1000 * ld, 1000 * radius, na2, s=1)

plt.show()
