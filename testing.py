import propagators as prop
import numpy as np
import matplotlib.pyplot as plt
import constants as cns

ra = 39700e3    # cns.re + 35000e3
rp = 600e3      # cns.re + 200e3

a = (ra + rp)/2
e = (ra - rp)/(ra + rp)
i = 63.4
lan = 45
w = 0
ta = 0

T = np.sqrt(4*np.pi**2*a**3/cns.mu)
x, y, z = prop.kepler_propagation(a, e, i, lan, w, ta, T*3, n=10000, j2=True)

ax = plt.axes(projection='3d')
ax.plot(x/1000, y/1000)
# plt.axis('equal')
ax.set_xlim([-a, a])
ax.set_ylim([-a, a])
ax.set_zlim([-a, a])
plt.show()
