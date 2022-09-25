import resonator as rs
import numpy as np
import matplotlib.pyplot as plt


NUM_OF_LENGTH_STEPS = 100
NUM_OF_RADIUS_STEPS = 100
LENGTH_STEP = 0.00005
RADIUS_STEP = 0.01

ld = np.zeros(NUM_OF_LENGTH_STEPS * NUM_OF_RADIUS_STEPS)
radius = np.zeros(NUM_OF_LENGTH_STEPS * NUM_OF_RADIUS_STEPS)
na1 = np.zeros(NUM_OF_LENGTH_STEPS * NUM_OF_RADIUS_STEPS)
na2 = np.zeros(NUM_OF_LENGTH_STEPS * NUM_OF_RADIUS_STEPS)
for i in range(NUM_OF_LENGTH_STEPS):
    for j in range(NUM_OF_RADIUS_STEPS):
        indx = i * NUM_OF_RADIUS_STEPS + j
        ld[indx] = (i - np.around(NUM_OF_LENGTH_STEPS / 2)) * LENGTH_STEP
        radius[indx] = (j - np.around(NUM_OF_RADIUS_STEPS / 2)) * RADIUS_STEP
        tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 0.005), rs.Mirror(0.008, 15, 0.010), rs.Mirror(0.2585 + ld[indx], 0, 0.250 + radius[indx]))
        na1[indx] = tmp_res.waist_search(1064 * 10 ** -9)[0, 0, 2]
        na2[indx] = tmp_res.waist_search(1064 * 10 ** -9)[1, 0, 2]
fig = plt.figure()
ax = plt.axes(projection= '3d')
ax.set_title("Numerical aperture (NA)")
ax.scatter3D(ld, radius, na1, s=1)

plt.show()
