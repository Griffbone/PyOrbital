import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import dev_functions
import functions as func
import dev_functions as dev
from scipy.optimize import minimize
import propagators

def yline(val, max_x):
    plt.plot([0, max_x], [val, val], 'k--')


def xline(val, max_y):
    plt.plot([val, val], [0, max_y], 'k--')


def circle(r, cx, cy, style):
    ts = np.linspace(0, 2*np.pi, 100)
    xs = r*np.cos(ts) + cx
    ys = r*np.sin(ts) + cy

    plt.plot(xs, ys, style)


ra = 0.3744e9
rp = 200000 + cns.re
soi = 66.1e6

sat = func.Elements((ra + rp)/2, (ra - rp)/(ra + rp), 0, 0, 0, 0, cns.mu)
moon = func.Elements(0.3844e9, 0, 0, 0, 0, 140, cns.mu)

#############################################################################
ts, xm, ym, zm = propagators.kepler_propagation(moon.a, moon.e, moon.i, 0, 0, moon.ta, sat.T/2)
_, xs, ys, zs = propagators.kepler_propagation(sat.a, sat.e, sat.i, 0, 0, sat.ta, max(ts))
rs = np.sqrt((xs - xm)**2 + (ys - ym)**2 + (zs - zm)**2)

tas, tam, _ = (dev_functions.patch_point(sat, moon, 66.1e6))
Dt = func.t_between(sat.a, sat.e, sat.ta, tas, sat.mu)

print(tas)

x, y, _ = sat.get_eci(tas)
plt.scatter(x, y)

x, y, _ = moon.get_eci(tam)
plt.scatter(x, y)
circle(soi, x, y, 'k--')

sat.perifocal_plot()
moon.perifocal_plot()

plt.show()
