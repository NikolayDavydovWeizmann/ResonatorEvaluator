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
L_3 = np.arange(Decimal('2'), Decimal('1000'), Decimal('2'))

f = (Decimal('1') + z_R ** Decimal('2')) / (Decimal('1') + ((Decimal('1') + z_R ** Decimal('2')) / L_3) + np.sqrt(((Decimal('1') + z_R ** Decimal('2')) / L_3)) ** Decimal('2') - z_R ** Decimal('2'))
z_R_1 = z_R * (Decimal('1') + z_R ** Decimal('2')) / (((Decimal('1') + z_R ** Decimal('2')) / f - Decimal('1')) ** Decimal('2') + z_R ** Decimal('2'))

r_1 = l_1 + z_R ** Decimal('2') / l_1
r_2 = Decimal('2') * f
r_3 = L_3 / Decimal('2') + z_R_1 ** Decimal('2') * Decimal('2') / L_3

res1 = rs.Resonator(3, rs.Mirror(Decimal('0'), Decimal('0'), r_1), rs.Mirror(Decimal('2') * l_1, Decimal('0'), r_2[0]), rs.Mirror(Decimal('2') * l_1 + L_3[0], Decimal('0'), r_3[0]))
print(res1.waist_search(FUND_WAVE))
na = np.zeros(L_3.size)
for i in range(L_3.size):
    res1 = rs.Resonator(3, rs.Mirror(Decimal('0'), Decimal('0'), r_1, in_plane_angle_deviation= C.power(Decimal('10'), Decimal('-13'))), rs.Mirror(Decimal('2') * l_1, Decimal('0'), r_2[i]), rs.Mirror(Decimal('2') * l_1 + L_3[i], Decimal('0'), r_3[i]))
    na[i] = res1.realign()

res2 = rs.Resonator(2, rs.Mirror(Decimal('0'), Decimal('0'), r_1, in_plane_angle_deviation= C.power(Decimal('10'), Decimal('-13'))), rs.Mirror(Decimal('2') * l_1, Decimal('0'), r_1))
na2 = res2.realign()
na = na2 / na
plt.plot(L_3 / 2 / l_1, na)
plt.show()
