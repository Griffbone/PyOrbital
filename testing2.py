import numpy as np
import functions as func
import constants as cns
import matplotlib.pyplot as plt


e = 2.7696
nu = np.deg2rad(100)

F = np.arctanh(np.sqrt((e - 1)/(e + 1))*np.tan(nu/2))*2
M = e*np.sinh(F) - F
t = ((h**3/mu_m**2)/(e**2 - 1)**(3/2))*M
print(F)

# # ra = 21000000
# # rp = 9600000
# ra = 384000e3
# rp = 200000 + cns.re
#
# a = (ra + rp)/2
# e = (ra - rp)/(ra + rp)
#
# sat = func.Elements(a, e, 0, 0, 0, 0, cns.mu)
#
# print(func.t_between(sat.a, sat.e, 0, 120, sat.mu))
# print(func.ta_change(sat.a, sat.e, 0, 3*60*60, sat.mu))
#
# ts = np.linspace(0, 3*sat.T, 100)
# tas = []
# for t in ts:
#     tas.append(func.ta_change(sat.a, sat.e, 0, t, sat.mu))
#
# plt.plot(ts, tas)
# plt.show()
