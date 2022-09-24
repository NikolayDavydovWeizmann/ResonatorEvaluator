import resonator as rs
import numpy as np
import matplotlib.pyplot as plt

m1 = rs.Mirror(0, 0, 0.005)
m2 = rs.Mirror(0.008, 30, 0.010)
m3 = rs.Mirror(0.2535, 0, 0.250)

res1 = rs.Resonator(3, m1, m2, m3)
print(res1.is_g_stable_sagittal())
print(res1.is_g_stable_tangential())
print(res1.waist_search(1064 * 10 ** -9))

m4 = rs.Mirror(0, 0, 0.005)
m5 = rs.Mirror(0.0091, 0, 0.005)
res2 = rs.Resonator(2, m4, m5)
print(res2.is_g_stable_sagittal())
print(res2.is_g_stable_tangential())
print(res2.waist_search(1064 * 10 ** -9))

tmp_res = rs.Resonator(2, rs.Mirror(0, 0, 0.005), rs.Mirror(0.0091, 0, 0.005))
init_split1 = tmp_res.longitude_split()
init_split2 = tmp_res.transverse_split()
print(init_split1)
print(init_split2)
angle = np.zeros(1000)
stblty1 = np.zeros(1000)
split1 = np.zeros(1000)
split2 = np.zeros(1000)
for i in range(1000):
    angle[i] = i * 0.001
    tmp_res = rs.Resonator(2, rs.Mirror(0, angle[i], 0.005), rs.Mirror(0.0091, 0, 0.005))
    tmp_res.realign()
    stblty1[i] = tmp_res.is_g_stable_sagittal()
    split1[i] = tmp_res.longitude_split() / init_split1 - 1
    split2[i] = tmp_res.transverse_split() / init_split2 - 1
plt.figure(1)
plt.plot(angle, stblty1)
plt.figure(2)
plt.plot(angle, split1)
plt.figure(3)
plt.plot(angle, split2)

tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 0.005), rs.Mirror(0.008, 30, 0.010), rs.Mirror(0.2535, 0, 0.250))
init_split1 = tmp_res.longitude_split()
init_split2 = tmp_res.transverse_split()
print(init_split1)
print(init_split2)
stblty2 = np.zeros(1000)
for i in range(1000):
    angle[i] = i * 0.001
    tmp_res = rs.Resonator(3, rs.Mirror(0, angle[i], 0.005), rs.Mirror(0.008, 30, 0.010), rs.Mirror(0.2535, 0, 0.250))
    tmp_res.realign()
    stblty1[i] = tmp_res.is_g_stable_sagittal()
    stblty2[i] = tmp_res.is_g_stable_tangential()
    split1[i] = tmp_res.longitude_split() / init_split1 - 1
    split2[i] = tmp_res.transverse_split() / init_split2 - 1
plt.figure(4)
plt.plot(angle, stblty1)
plt.plot(angle, stblty2)
plt.figure(5)
plt.plot(angle, split1)
plt.figure(6)
plt.plot(angle, split2)

plt.show()
