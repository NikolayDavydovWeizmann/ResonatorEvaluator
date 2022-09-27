import resonator as rs
import numpy as np
import matplotlib.pyplot as plt

ANGLE_STEP = 4 * 180 / np.pi / 1000000 / 1000
NUMBER_OF_STEPS = 1000

tmp_res = rs.Resonator(2, rs.Mirror(0, 0, 0.005), rs.Mirror(0.00999956, 0, 0.005))
init_split1 = tmp_res.longitude_split()
init_split2 = tmp_res.transverse_split()
fund_wave = tmp_res.fund_lambda_choice(1064 * 10 ** -9)
init_na1 = tmp_res.waist_search(fund_wave)[0, 0, 2]
print(init_split1)
print(init_split2)
print(fund_wave)
angle = np.zeros(NUMBER_OF_STEPS)
stblty1 = np.zeros(NUMBER_OF_STEPS)
split1 = np.zeros(NUMBER_OF_STEPS)
split2 = np.zeros(NUMBER_OF_STEPS)
wave_length = np.zeros(NUMBER_OF_STEPS)
na1 = np.zeros(NUMBER_OF_STEPS)
rot = np.zeros(NUMBER_OF_STEPS)
for i in range(NUMBER_OF_STEPS):
    angle[i] = (i - np.around(NUMBER_OF_STEPS / 2)) * ANGLE_STEP
    tmp_res = rs.Resonator(2, rs.Mirror(0, angle[i], 0.005), rs.Mirror(0.00999956, 0, 0.005))
    rot[i] = tmp_res.realign()
    wave_length[i] = tmp_res.fund_lambda_choice(1064 * 10 ** -9) / fund_wave - 1
    stblty1[i] = tmp_res.is_g_stable_sagittal()
    split1[i] = tmp_res.longitude_split() / init_split1 - 1
    split2[i] = tmp_res.transverse_split() / init_split2 - 1
    na1[i] = tmp_res.waist_search(fund_wave)[0, 0, 2] / init_na1 - 1
plt.figure(1)
plt.title("Numerical aperture (NA)")
plt.xlabel("Angle deviation")
plt.ylabel("Relative deviation of NA")
plt.plot(angle, na1)
"""
plt.figure(2)
plt.title("Longitude split (LS)")
plt.xlabel("Angle deviation")
plt.ylabel("Relative deviation of LS")
plt.plot(angle, split1)
plt.figure(3)
plt.title("Transverse split (TS)")
plt.xlabel("Angle deviation")
plt.ylabel("Relative deviation of TS")
plt.plot(angle, split2)
"""
plt.figure(4)
plt.title("Fundamental wave length (FWL)")
plt.xlabel("Angle deviation")
plt.ylabel("Relative deviation of FWL")
plt.plot(angle, wave_length)
plt.figure(5)
plt.title("Rotation of axis of a mode")
plt.xlabel("Angle deviation")
plt.ylabel("Angle of rotation")
plt.plot(angle, rot)


tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 0.005), rs.Mirror(0.0099999, 10, 0.010), rs.Mirror(0.300, 0, 0.595))
init_split1 = tmp_res.longitude_split()
init_split2 = tmp_res.transverse_split()
fund_wave = tmp_res.fund_lambda_choice(1064 * 10 ** -9)
init_na1 = tmp_res.waist_search(fund_wave)[0, 0, 2]
init_na2 = tmp_res.waist_search(fund_wave)[0, 1, 2]
print(init_split1)
print(init_split2)
print(fund_wave)
stblty2 = np.zeros(NUMBER_OF_STEPS)
na2 = np.zeros(NUMBER_OF_STEPS)
for i in range(NUMBER_OF_STEPS):
    angle[i] = (i - np.around(NUMBER_OF_STEPS / 2)) * ANGLE_STEP
    tmp_res = rs.Resonator(3, rs.Mirror(0, angle[i], 0.005), rs.Mirror(0.0099999, 10, 0.010), rs.Mirror(0.300, 0, 0.595))
    rot[i] = tmp_res.realign()
    wave_length[i] = tmp_res.fund_lambda_choice(1064 * 10 ** -9) / fund_wave - 1
    stblty1[i] = tmp_res.is_g_stable_sagittal()
    stblty2[i] = tmp_res.is_g_stable_tangential()
    split1[i] = tmp_res.longitude_split() / init_split1 - 1
    split2[i] = tmp_res.transverse_split() / init_split2 - 1
    na1[i] = tmp_res.waist_search(fund_wave)[0, 0, 2] / init_na1 - 1
    na2[i] = tmp_res.waist_search(fund_wave)[0, 1, 2] / init_na2 - 1
plt.figure(6)
plt.title("Numerical aperture (NA)")
plt.xlabel("Angle deviation")
plt.ylabel("Relative deviation of NA")
plt.plot(angle, na1)
plt.plot(angle, na2)
"""
plt.figure(7)
plt.title("Longitude split (LS)")
plt.xlabel("Angle deviation")
plt.ylabel("Relative deviation of LS")
plt.plot(angle, split1)
plt.figure(8)
plt.title("Transverse split (TS)")
plt.xlabel("Angle deviation")
plt.ylabel("Relative deviation of TS")
plt.plot(angle, split2)
"""
plt.figure(9)
plt.title("Fundamental wave length (FWL)")
plt.xlabel("Angle deviation")
plt.ylabel("Relative deviation of FWL")
plt.plot(angle, wave_length)
plt.figure(10)
plt.title("Rotation of axis of a mode")
plt.xlabel("Angle deviation")
plt.ylabel("Angle of rotation")
plt.plot(angle, rot)

for i in range(NUMBER_OF_STEPS):
    angle[i] = (i - np.around(NUMBER_OF_STEPS / 2)) * ANGLE_STEP
    tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 0.005), rs.Mirror(0.0099999, 10, 0.010), rs.Mirror(0.300, angle[i], 0.595))
    rot[i] = tmp_res.realign()
    wave_length[i] = tmp_res.fund_lambda_choice(1064 * 10 ** -9) / fund_wave - 1
    stblty1[i] = tmp_res.is_g_stable_sagittal()
    stblty2[i] = tmp_res.is_g_stable_tangential()
    split1[i] = tmp_res.longitude_split() / init_split1 - 1
    split2[i] = tmp_res.transverse_split() / init_split2 - 1
    na1[i] = tmp_res.waist_search(fund_wave)[0, 0, 2] / init_na1 - 1
    na2[i] = tmp_res.waist_search(fund_wave)[0, 1, 2] / init_na2 - 1
plt.figure(11)
plt.title("Numerical aperture (NA)")
plt.xlabel("Angle deviation")
plt.ylabel("Relative deviation of NA")
plt.plot(angle, na1)
plt.plot(angle, na2)
"""
plt.figure(12)
plt.title("Longitude split (LS)")
plt.xlabel("Angle deviation")
plt.ylabel("Relative deviation of LS")
plt.plot(angle, split1)
plt.figure(13)
plt.title("Transverse split (TS)")
plt.xlabel("Angle deviation")
plt.ylabel("Relative deviation of TS")
plt.plot(angle, split2)
"""
plt.figure(14)
plt.title("Fundamental wave length (FWL)")
plt.xlabel("Angle deviation")
plt.ylabel("Relative deviation of FWL") 
plt.plot(angle, wave_length)
plt.figure(15)
plt.title("Rotation of axis of a mode")
plt.xlabel("Angle deviation")
plt.ylabel("Angle of rotation")
plt.plot(angle, rot)

plt.show()
