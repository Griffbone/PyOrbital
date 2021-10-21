import propagators as prop
import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import cartopy.crs as crs
import astrotime as at
import functions as func
from tle import read_tle


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


# Get orbital elements from the TLE
l1 = 'ISS (ZARYA)'
l2 = '1 25544U 98067A   21292.89044074  .00004617  00000-0  92733-4 0  9998'
l3 = '2 25544  51.6431  86.1166 0004166 125.7553 304.8177 15.48743137307929'
iss = read_tle(l1, l2, l3)

epoch = iss['epoch']
a = iss['sma']
e = iss['eccentricity']
i = iss['inclination']
lan = iss['raan']
w = iss['argpe']
ma = iss['ma']
ta = func.kepler_ta(e, ma)

# get position of ISS at epoch
x, y, z = func.elements_to_vector(a, e, i, lan, w, ta)
print(np.array([x, y, z])/1000)

t, x, y, z = prop.kepler_propagation(a, e, i, lan, w, ta, 60*60*5, n=1000)
# phi, lam = func.eci_to_ll(x, y, z, jdn)
#
print(np.array([x[0], y[0], z[0]])/1000)
# print(np.array([x[-1], y[-1], z[-1]])/1000)
