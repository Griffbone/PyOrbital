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


# l1 = 'ISS (ZARYA)'
# l2 = '1 25544U 98067A   21294.13453560  .00002095  00000-0  46566-4 0  9992'
# l3 = '2 25544  51.6430  79.9628 0004254 128.2192  43.3615 15.48750324308'
# iss = read_tle(l1, l2, l3)
#
# epoch = iss['epoch']
# a = iss['sma']
# e = iss['eccentricity']
# i = iss['inclination']
# lan = iss['raan']
# w = iss['argpe']
# ma = iss['ma']
# ta = func.kepler_ta(e, ma)
#
#
# j0 = at.datetime_to_jd(at.tle_to_datetime(21294.13453560))
# now = at.datetime_to_jd(datetime.utcnow())
# dt = (now - j0) * 60**2 * 24
#
# t, x, y, z = prop.kepler_propagation(a, e, i, lan, w, ta, dt, n=1000, j2=True)
# t2, x2, y2, z2 = prop.kepler_propagation(a, e, i, lan, w, ta, dt, n=1000, j2=False)
#
# lats = []
# lons = []
# lats2 = []
# lons2 = []
#
# for i in range(0, len(t)):
#     jdn = j0 + t[i]/(60**2 * 24)
#     lat, lon = func.eci_to_ll(x[i], y[i], z[i], jdn)
#
#     lats.append(lat)
#     lons.append(lon)
#
#     lat, lon = func.eci_to_ll(x2[i], y2[i], z2[i], jdn)
#     lats2.append(lat)
#     lons2.append(lon)
#
# ax = plt.axes(projection=crs.PlateCarree())
# ax.coastlines(resolution='110m')
#
# ax.scatter(lons[-1], lats[-1], marker='x', color='r', s=100)
# # ax.plot(lons, lats, transform=crs.Geodetic())
#
# ax.scatter(lons2[-1], lats2[-1], marker='x', color='b', s=100)
# # ax.plot(lons2, lats2, transform=crs.Geodetic())
#
# ax.set_ylim([-180, 180])
#
# ax.add_feature(Nightshade(datetime.utcnow(), alpha=0.2))
# plt.show()
