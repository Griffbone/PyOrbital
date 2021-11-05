import numpy as np
import functions as func
import constants as cns
import matplotlib.pyplot as plt

# ra = 21000000
# rp = 9600000
ra = 384000e3
rp = 200000 + cns.re

a = (ra + rp)/2
e = (ra - rp)/(ra + rp)

sat = func.Elements(a, e, 0, 0, 0, 0, cns.mu)

print(func.t_between(sat.a, sat.e, 0, 120, sat.mu))
print(func.ta_change(sat.a, sat.e, 0, 3*60*60, sat.mu))

ts = np.linspace(0, 3*sat.T, 100)
tas = []
for t in ts:
    tas.append(func.ta_change(sat.a, sat.e, 0, t, sat.mu))

plt.plot(ts, tas)
plt.show()
