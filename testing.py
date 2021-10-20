import propagators as prop
import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import cartopy.crs as crs

#
import astrotime as at
import functions as func


# ISS (ZARYA)
# 1 25544U 98067A   21292.89044074  .00004617  00000-0  92733-4 0  9998
# 2 25544  51.6431  86.1166 0004166 125.7553 304.8177 15.48743137307929


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


a = 6798.515912586929e3
e = 0.0004166
i = 51.6431
lan = 86.1166
w = 125.7553
ta = func.kepler_ta(e, 304.8177)

t0 = at.date_to_jd(19, 10, 2021) + 21/24 + 22/(60*24) + 14/(60*60*24)
t = at.date_to_jd(20, 10, 2021) + 14/24 + 8/(60*24) + 22/(60**2*24)

dt = (t-t0)*24
hr = int(dt)
mi = int((dt - hr)*60)
sec = (dt*60**2 - hr*60**2 - mi*60)
print('Time since epoch: {}:{}:{:.2f}'.format(hr, mi, sec))

ts, x, y, z = prop.kepler_propagation(a, e, i, lan, w, ta, (t-t0)*60**2*24, n=1000, j2=False)

print(np.array([x[-1], y[-1], z[-1]])/1000)

# ts = ts + (t - t0)*24*60**2
# lat, lon = ground_track(at.date_to_jd(19, 10, 2021), ts, x, y, z)
#
# ax = plt.axes(projection=crs.PlateCarree())
# ax.coastlines()
# ax.plot(lon, lat, transform=crs.Geodetic())
#
# ax.scatter(lon[-1], lat[-1])
#
#
# ts, x, y, z = prop.kepler_propagation(a, e, i, lan, w, ta, (t-t0)*60**2*24, n=1000, j2=False)
# ts = ts + (t - t0)*24*60**2
# lat, lon = ground_track(at.date_to_jd(19, 10, 2021), ts, x, y, z)
#
#
# ax.plot(lon, lat, transform=crs.Geodetic())
# ax.scatter(lon[-1], lat[-1])
#
#
# plt.show()
