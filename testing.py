import propagators as prop
import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import cartopy.crs as crs

#
import astrotime as at
import functions as func


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

    for i in range(0, len(xs)):
        x = xs[i]
        y = ys[i]
        z = zs[i]

        r = np.sqrt(x**2 + y**2 + z**2)
        phi = np.arcsin(z/r)
        lam = np.arctan2(y, x)

        lam = np.degrees(lam)
        tg = at.theta_g(t0, ts[i]/3600)

        lam = lam + tg - int((lam+tg)/360)*360

        phis.append(np.degrees(phi))
        lams.append(lam)

    return phis, lams


a = (cns.re*2 + 417e3 + 423e3)/2
e = 0.0004158
i = 51.6433
lan = 88.0366
w = 127.3682
ta = func.kepler_ta(e, 299.1142)

T = np.sqrt(4*np.pi**2*a**3/cns.mu)
t0 = at.date_to_jd(19, 10, 2021) + 12/24 + 3/(60*24) + 40/(60*60*24)
t = at.date_to_jd(19, 10, 2021) + (9+12)/12 + 52/(60*24)

ts, x, y, z = prop.kepler_propagation(a, e, i, lan, w, ta, (t-t0)*60**2*24, n=1000, j2=False)

lat, lon = ground_track(at.date_to_jd(19, 10, 2021), ts, x, y, z)
# plt.plot(lat)
# plt.plot(lon)
# plt.show()

ax = plt.axes(projection=crs.PlateCarree())
ax.coastlines()

ax.plot(lon, lat, transform=crs.Geodetic())
ax.scatter(lon[-1], lat[-1])
plt.show()
