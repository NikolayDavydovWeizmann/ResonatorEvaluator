from decimal import *
import resonator as rs
import numpy as np
import matplotlib.pyplot as plt

C = getcontext()
C.prec = 28
PI = Decimal('3.141592653589793238462643383279')

FUND_WAVE = Decimal('1064') * C.power(Decimal('10'), Decimal('-7'))

"""
ANGLE = Decimal('15') / Decimal('180') * PI
"""

NA = Decimal('0.2')
z_R = FUND_WAVE / PI / C.power(NA, Decimal('2'))

l_1 = Decimal('1')
L_3 = np.arange(Decimal('2'), Decimal('200'), Decimal('2'))

f = np.empty(L_3.size, dtype='O')
for i in range(L_3.size):
    f[i] = (Decimal('1') + z_R ** Decimal('2')) / (Decimal('1') + ((Decimal('1') + z_R ** Decimal('2')) / L_3[i]) + C.sqrt(((Decimal('1') + z_R ** Decimal('2')) / L_3[i]) ** Decimal('2') - z_R ** Decimal('2')))
z_R_1 = z_R * (Decimal('1') + z_R ** Decimal('2')) / (((Decimal('1') + z_R ** Decimal('2')) / f - Decimal('1')) ** Decimal('2') + z_R ** Decimal('2'))

r_1 = l_1 + z_R ** Decimal('2') / l_1
r_2 = Decimal('2') * f
r_3 = L_3 / Decimal('2') + z_R_1 ** Decimal('2') * Decimal('2') / L_3

na = np.empty(L_3.size, dtype= 'O')
for i in range(L_3.size):
    res1 = rs.Resonator(3, rs.Mirror(Decimal('0'), Decimal('0'), r_1), rs.Mirror(Decimal('2') * l_1, Decimal('0'), r_2[i]), rs.Mirror(Decimal('2') * l_1 + L_3[i], Decimal('0'), r_3[i]))
    na_tmp = res1.waist_search(FUND_WAVE)[0, 0, 2]
    res1 = rs.Resonator(3, rs.Mirror(Decimal('0'), Decimal('0'), r_1), rs.Mirror(Decimal('2') * l_1 + Decimal('10') ** Decimal('-10'), Decimal('0'), r_2[i]), rs.Mirror(Decimal('2') * l_1 + Decimal('10') ** Decimal('-10') + L_3[i], Decimal('0'), r_3[i]))
    na[i] = res1.waist_search(FUND_WAVE)[0, 0, 2] / na_tmp - Decimal('1')

res2 = rs.Resonator(2, rs.Mirror(Decimal('0'), Decimal('0'), r_1), rs.Mirror(Decimal('2') * l_1, Decimal('0'), r_1))
na2_init = res2.waist_search(FUND_WAVE)[0, 0, 2]
res2 = rs.Resonator(2, rs.Mirror(Decimal('0'), Decimal('0'), r_1), rs.Mirror(Decimal('2') * l_1 + Decimal('10') ** Decimal('-10'), Decimal('0'), r_1))
na2 = res2.waist_search(FUND_WAVE)[0, 0, 2] / na2_init - Decimal('1')
plt.plot(L_3 / Decimal('2') / l_1, na / na2)
plt.yscale('log')
plt.show()
