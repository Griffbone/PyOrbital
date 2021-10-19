import propagators as prop
import constants as cns
import numpy as np
import matplotlib.pyplot as plt

ap = cns.re + 35000000
pe = cns.re + 200000

a = (ap + pe)/2
e = (ap - pe)/(ap + pe)

i = 28
lan = 0
w = 0
ta = 0

T = np.sqrt(4*np.pi**2*a**3/cns.mu)
ts, tas = prop.kepler_propagation(a, e, i, lan, w, ta, T)

plt.plot(ts, tas*180/np.pi)
plt.show()
