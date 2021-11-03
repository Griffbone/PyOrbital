import propagators as prop
import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import cartopy.crs as crs
import astrotime as at
import functions as func
from tle import read_tle
from cartopy.feature.nightshade import Nightshade
from datetime import datetime, timezone, timedelta


def ground_track(t0, ts, xs, ys, zs):
    """ Function to return geodetic lat/lon from ECI positions
        :param t0: epoch (JDN)
        :param ts: array of times since epoch (s)
        :param xs: ECI X position
        :param ys: ECI Y position
        :param zs: ECI Z position
    """

    phis = []
    lams = []

    for i in range(0, len(ts)):
        t = ts[i]/60**2
        phi, lam = func.eci_to_ll(xs[i], ys[i], zs[i], t0, t)

        phis.append(phi)
        lams.append(lam)

    return phis, lams


rp = 200e3 + cns.re
ra = 485000e3 + cns.re
a = (rp + ra)/2
e = (ra - rp)/(ra + rp)

T = 1000000

ts, x1, y1, z1 = prop.kepler_propagation(a, e, 0, 0, 0, 0, T, n=1000, j2=False)

ts, x2, y2, z2 = prop.kepler_propagation(ra, 0, 0, 0, 0, 90, T, n=1000, j2=False)

dx = x2 - x1
dy = y2 - y1
dz = z2 - z1

rs = np.sqrt(dx**2 + dy**2 + dz**2)

# for i in range(0, len(ts)):
#     r = np.array([x1])
#     r = np.sqrt(x[])

# plt.axis('equal')
plt.plot(ts, rs)
plt.show()
