from decimal import *
import numpy as np
import matplotlib.pyplot as plt

C = getcontext()
C.prec = 28
PI = Decimal('3.141592653589793238462643383279')

FUND_WAVE = Decimal('1064') * C.power(Decimal('10'), Decimal('-7'))

NA = Decimal('0.2')
z_R = FUND_WAVE / PI / C.power(NA, Decimal('2'))

l_1 = Decimal('1')
L_3 = np.arange(Decimal('2'), Decimal('1000'), Decimal('2'))

f = np.empty(L_3.size, dtype='O')
for i in range(L_3.size):
    f[i] = (Decimal('1') + z_R ** Decimal('2')) / (Decimal('1') + ((Decimal('1') + z_R ** Decimal('2')) / L_3[i]) + C.sqrt(((Decimal('1') + z_R ** Decimal('2')) / L_3[i]) ** Decimal('2') - z_R ** Decimal('2')))

z_R_1 = z_R * (Decimal('1') + z_R ** Decimal('2')) / (((Decimal('1') + z_R ** Decimal('2')) / f - Decimal('1')) ** Decimal('2') + z_R ** Decimal('2'))
l_0_1 = z_R_1 / z_R * ((Decimal('1') + z_R ** Decimal('2')) / f - Decimal('1'))


xi_1 = Decimal('1') - z_R ** Decimal('2') / l_1
xi_2 = l_0_1 - Decimal('2') * z_R_1 ** Decimal('2') / L_3
r_1 = l_1 + z_R ** Decimal('2') / l_1
r_3 = L_3 / Decimal('2') + z_R_1 ** Decimal('2') * Decimal('2') / L_3

acc = (f * xi_2 - xi_1 * (xi_2 - f)) / (r_1 * (xi_2 - f))
acc3 = (f * xi_2 - xi_1 * (xi_2 - f)) / r_3 / f


acc_2mirror = Decimal('2') * (Decimal('1') - Decimal('1') / r_1)

plt.figure(1)
plt.plot(L_3 / Decimal('2'), acc / acc_2mirror)
plt.figure(2)
plt.plot(L_3 / Decimal('2'), acc3 / acc_2mirror)
plt.show()
