import resonator as rs
import numpy as np
import matplotlib.pyplot as plt

tmp_res = rs.Resonator(2, rs.Mirror(0, 0, 0.005), rs.Mirror(0.0091, 0, 0.005))
init_split1 = tmp_res.longitude_split()
init_split2 = tmp_res.transverse_split()
fund_wave = tmp_res.fund_lambda_choice(1064 * 10 ** -9)
print(init_split1)
print(init_split2)
print(fund_wave)
angle = np.zeros(1000)
stblty1 = np.zeros(1000)
split1 = np.zeros(1000)
split2 = np.zeros(1000)
wave_length = np.zeros(1000)
for i in range(1000):
    angle[i] = i * 0.0001
    tmp_res = rs.Resonator(2, rs.Mirror(0, angle[i], 0.005), rs.Mirror(0.0091, 0, 0.005))
    tmp_res.realign()
    wave_length[i] = tmp_res.fund_lambda_choice(1064 * 10 ** -9) / fund_wave - 1
    stblty1[i] = tmp_res.is_g_stable_sagittal()
    split1[i] = tmp_res.longitude_split() / init_split1 - 1
    split2[i] = tmp_res.transverse_split() / init_split2 - 1
plt.figure(2)
plt.plot(angle, split1)
plt.figure(3)
plt.plot(angle, split2)
plt.figure(4)
plt.plot(angle, wave_length)


tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 0.005), rs.Mirror(0.008, 30, 0.010), rs.Mirror(0.2535, 0, 0.250))
init_split1 = tmp_res.longitude_split()
init_split2 = tmp_res.transverse_split()
fund_wave = tmp_res.fund_lambda_choice(1064 * 10 ** -9)
print(init_split1)
print(init_split2)
print(fund_wave)
stblty2 = np.zeros(1000)
for i in range(1000):
    angle[i] = i * 0.0001
    tmp_res = rs.Resonator(3, rs.Mirror(0, angle[i], 0.005), rs.Mirror(0.008, 30, 0.010), rs.Mirror(0.2535, 0, 0.250))
    tmp_res.realign()
    wave_length[i] = tmp_res.fund_lambda_choice(1064 * 10 ** -9) / fund_wave - 1
    stblty1[i] = tmp_res.is_g_stable_sagittal()
    stblty2[i] = tmp_res.is_g_stable_tangential()
    split1[i] = tmp_res.longitude_split() / init_split1 - 1
    split2[i] = tmp_res.transverse_split() / init_split2 - 1
plt.figure(6)
plt.plot(angle, split1)
plt.figure(7)
plt.plot(angle, split2)
plt.figure(8)
plt.plot(angle, wave_length)

plt.show()
