from math import nan
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt
 
plt.rcParams.update({'font.size': 30})

NUMBER_OF_STEPS = 10000
FUND_WAVE = 1064 * 10 ** -7

na1 = np.zeros(NUMBER_OF_STEPS)
na2 = np.zeros(NUMBER_OF_STEPS)
chi = np.zeros(NUMBER_OF_STEPS)

ANGLE = 0 / 180 * np.pi
R_2 = 10 / 11
R_3 = 200 / 11
xi = 1.01

MIN_CHI = np.max(np.array([R_2 / 2, xi - np.cos(ANGLE), xi - 1 + R_2 / 2 / np.cos(ANGLE)])) * 0.9999999999999
MAX_CHI = np.cos(ANGLE) * 0.9999999999999
CHI_STEPS = (MAX_CHI - MIN_CHI) / (NUMBER_OF_STEPS - 1)

for i in range(NUMBER_OF_STEPS):
    chi[i] = MIN_CHI + i * CHI_STEPS
    tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 1 - R_2 / 2 / chi[i]), rs.Mirror(1, ANGLE, R_2), rs.Mirror(1 + R_3 + R_2 / 2 / (xi - chi[i]), 0, R_3))
    if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
        tmp_na1 = tmp_res.waist_search(FUND_WAVE)[0, 0, 2]
        tmp_na2 = tmp_res.waist_search(FUND_WAVE)[0, 1, 2]
        na1[i] = tmp_na1 if tmp_na1 < 0.21 else nan
        na2[i] = tmp_na2 if tmp_na2 < 0.21 else nan
    else:
        na1[i] = nan
        na2[i] = nan

plt.figure(1)
plt.plot(chi, na1, color= '#1f77b4', label= r'$R_2 = 0.91, R_3 = 18.18$')
plt.plot(chi, na2, color= '#1f77b4', linestyle= 'dashed')


R_2 = 8 / 11
MIN_CHI = np.max(np.array([R_2 / 2, xi - np.cos(ANGLE), xi - 1 + R_2 / 2 / np.cos(ANGLE)])) * 0.9999999999999
MAX_CHI = np.cos(ANGLE) * 0.9999999999999
CHI_STEPS = (MAX_CHI - MIN_CHI) / (NUMBER_OF_STEPS - 1)

for i in range(NUMBER_OF_STEPS):
    chi[i] = MIN_CHI + i * CHI_STEPS
    tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 1 - R_2 / 2 / chi[i]), rs.Mirror(1, ANGLE, R_2), rs.Mirror(1 + R_3 + R_2 / 2 / (xi - chi[i]), 0, R_3))
    if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
        tmp_na1 = tmp_res.waist_search(FUND_WAVE)[0, 0, 2]
        tmp_na2 = tmp_res.waist_search(FUND_WAVE)[0, 1, 2]
        na1[i] = tmp_na1 if tmp_na1 < 0.21 else nan
        na2[i] = tmp_na2 if tmp_na2 < 0.21 else nan
    else:
        na1[i] = nan
        na2[i] = nan

plt.plot(chi, na1, color= '#ff7f0e', label= r'$R_2 = 0.73, R_3 = 18.18$')
plt.plot(chi, na2, color= '#ff7f0e', linestyle= 'dashed')


R_2 = 6 / 11
MIN_CHI = np.max(np.array([R_2 / 2, xi - np.cos(ANGLE), xi - 1 + R_2 / 2 / np.cos(ANGLE)])) * 0.9999999999999
MAX_CHI = np.cos(ANGLE) * 0.9999999999999
CHI_STEPS = (MAX_CHI - MIN_CHI) / (NUMBER_OF_STEPS - 1)

for i in range(NUMBER_OF_STEPS):
    chi[i] = MIN_CHI + i * CHI_STEPS
    tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 1 - R_2 / 2 / chi[i]), rs.Mirror(1, ANGLE, R_2), rs.Mirror(1 + R_3 + R_2 / 2 / (xi - chi[i]), 0, R_3))
    if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
        tmp_na1 = tmp_res.waist_search(FUND_WAVE)[0, 0, 2]
        tmp_na2 = tmp_res.waist_search(FUND_WAVE)[0, 1, 2]
        na1[i] = tmp_na1 if tmp_na1 < 0.21 else nan
        na2[i] = tmp_na2 if tmp_na2 < 0.21 else nan
    else:
        na1[i] = nan
        na2[i] = nan

plt.plot(chi, na1, color= '#2ca02c', label= r'$R_2 = 0.55, R_3 = 18.18$')
plt.plot(chi, na2, color= '#2ca02c', linestyle= 'dashed')

"""plt.yscale('log')"""
plt.xlabel(r'$\chi$')
plt.ylabel('NA')
plt.grid(visible= True, which= 'both')
plt.legend()
plt.get_current_fig_manager().window.state('zoomed')

R_2 = 10 / 11
R_3 = 100 / 11

MIN_CHI = np.max(np.array([R_2 / 2, xi - np.cos(ANGLE), xi - 1 + R_2 / 2 / np.cos(ANGLE)])) * 0.9999999999999
MAX_CHI = np.cos(ANGLE) * 0.9999999999999
CHI_STEPS = (MAX_CHI - MIN_CHI) / (NUMBER_OF_STEPS - 1)

for i in range(NUMBER_OF_STEPS):
    chi[i] = MIN_CHI + i * CHI_STEPS
    tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 1 - R_2 / 2 / chi[i]), rs.Mirror(1, ANGLE, R_2), rs.Mirror(1 + R_3 + R_2 / 2 / (xi - chi[i]), 0, R_3))
    if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
        tmp_na1 = tmp_res.waist_search(FUND_WAVE)[0, 0, 2]
        tmp_na2 = tmp_res.waist_search(FUND_WAVE)[0, 1, 2]
        na1[i] = tmp_na1 if tmp_na1 < 0.21 else nan
        na2[i] = tmp_na2 if tmp_na2 < 0.21 else nan
    else:
        na1[i] = nan
        na2[i] = nan

plt.figure(2)
plt.plot(chi, na1, color= '#1f77b4', label= r'$R_2 = 0.91, R_3 = 9.09$')
plt.plot(chi, na2, color= '#1f77b4', linestyle= 'dashed')

R_3 = 200 / 11

MIN_CHI = np.max(np.array([R_2 / 2, xi - np.cos(ANGLE), xi - 1 + R_2 / 2 / np.cos(ANGLE)])) * 0.9999999999999
MAX_CHI = np.cos(ANGLE) * 0.9999999999999
CHI_STEPS = (MAX_CHI - MIN_CHI) / (NUMBER_OF_STEPS - 1)

for i in range(NUMBER_OF_STEPS):
    chi[i] = MIN_CHI + i * CHI_STEPS
    tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 1 - R_2 / 2 / chi[i]), rs.Mirror(1, ANGLE, R_2), rs.Mirror(1 + R_3 + R_2 / 2 / (xi - chi[i]), 0, R_3))
    if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
        tmp_na1 = tmp_res.waist_search(FUND_WAVE)[0, 0, 2]
        tmp_na2 = tmp_res.waist_search(FUND_WAVE)[0, 1, 2]
        na1[i] = tmp_na1 if tmp_na1 < 0.21 else nan
        na2[i] = tmp_na2 if tmp_na2 < 0.21 else nan
    else:
        na1[i] = nan
        na2[i] = nan

plt.plot(chi, na1, color= '#ff7f0e', label= r'$R_2 = 0.91, R_3 = 18.18$')
plt.plot(chi, na2, color= '#ff7f0e', linestyle= 'dashed')

R_3 = 300 / 11

MIN_CHI = np.max(np.array([R_2 / 2, xi - np.cos(ANGLE), xi - 1 + R_2 / 2 / np.cos(ANGLE)])) * 0.9999999999999
MAX_CHI = np.cos(ANGLE) * 0.9999999999999
CHI_STEPS = (MAX_CHI - MIN_CHI) / (NUMBER_OF_STEPS - 1)

for i in range(NUMBER_OF_STEPS):
    chi[i] = MIN_CHI + i * CHI_STEPS
    tmp_res = rs.Resonator(3, rs.Mirror(0, 0, 1 - R_2 / 2 / chi[i]), rs.Mirror(1, ANGLE, R_2), rs.Mirror(1 + R_3 + R_2 / 2 / (xi - chi[i]), 0, R_3))
    if tmp_res.is_g_stable_sagittal() and tmp_res.is_g_stable_tangential():
        tmp_na1 = tmp_res.waist_search(FUND_WAVE)[0, 0, 2]
        tmp_na2 = tmp_res.waist_search(FUND_WAVE)[0, 1, 2]
        na1[i] = tmp_na1 if tmp_na1 < 0.21 else nan
        na2[i] = tmp_na2 if tmp_na2 < 0.21 else nan
    else:
        na1[i] = nan
        na2[i] = nan

plt.plot(chi, na1, color= '#2ca02c', label= r'$R_2 = 0.91, R_3 = 27.27$')
plt.plot(chi, na2, color= '#2ca02c', linestyle= 'dashed')

"""plt.yscale('log')"""
plt.xlabel(r'$\chi$')
plt.ylabel('NA')
plt.grid(visible= True, which= 'both')
plt.legend()
plt.get_current_fig_manager().window.state('zoomed')

plt.show()
