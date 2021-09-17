import matplotlib.pyplot as plt
import numpy as np
import pyOrbital.constants as cons
import pyOrbital.functions as func
import pyOrbital.astrotime as time
from scipy.optimize import fsolve
import time as timelib
import cartopy.crs as ccrs
import pyOrbital.constants as cons
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER


def draw_sphere(ax, r):
    u, v = np.mgrid[0:2 * np.pi:200j, 0:np.pi:100j]
    x = r * np.cos(u) * np.sin(v)
    y = r * np.sin(u) * np.sin(v)
    z = r * np.cos(v)

    ax.plot([0, 2*r], [0, 0], [0, 0], 'r')
    ax.plot([0, 0], [0, 2*r], [0, 0], 'g')
    ax.plot([0, 0], [0, 0], [0, 2*r], 'b')

    ax.plot_wireframe(x, y, z, color="k")


def elements_to_vec(a, e, i, omega, w, ta):
    """ Function to get ECI position vector from orbital elements
    """

    # secondary orbit geometry
    b = a*np.sqrt(1 - e**2)
    p = b**2/a

    # perifocal coordinates
    r = p/(1 + e*np.cos(ta))
    x = r*np.cos(ta)
    y = r*np.sin(ta)
    z = 0

    # rotation matrices (perifocal to ECI)
    rvec = np.array([x, y, z])

    Z1 = func.rotmat(omega, 'z')
    X1 = func.rotmat(i, 'x')
    Z2 = func.rotmat(w, 'z')

    R = np.matmul(Z1, X1)
    R = np.matmul(R, Z2)

    # ECI coordinates
    rvec_eci = np.matmul(R, rvec)

    return rvec_eci


def get_geo_coords(r, jd):
    """ Function to get geographic coordinates of a satellite from ECI position
        This function assumes a spherical Earth

        :param r: ECI position vector
        :param jd: julian data associated with instantaneous satellite position
    """

    # get right ascension and declination
    theta = np.arctan2(r[1], r[0])              # right ascension
    phi = np.arcsin(r[2]/np.linalg.norm(r))     # declination

    # get Greenwich sidereal time
    _, theta_g, _ = time.jd_to_t0(jd)
    theta_g = np.radians(theta_g)

    # calculate east longitude
    theta = theta - theta_g

    while theta < 0:
        theta += 2*np.pi

    return theta, phi


# INSPIRATION4
# 1 49220U 21084A   21259.48981751 -.00030321  00000-0 -25580-2 0  9997
# 2 49220  51.6453 248.7948 0006916 327.7800 175.2734 14.97540535    73


t0 = time.date_to_jd(16, 9, 2021, hour=11, mins=45, sec=20)
a = (659e3 + cons.re + 579e3 + cons.re)/2
e = .0006916
i = np.radians(51.6453)
omega = np.radians(249.8315)
w = np.radians(326.8812)
T = 2*np.pi*np.sqrt(a**3/cons.mu)

M = np.radians(175.2734)
t = (M/(2*np.pi))*T

t0 -= t/3600/24


jd = timelib.time() / 60 / 60 / 24 + time.date_to_jd(1, 1, 1970)
ts = np.linspace(t0, jd, 1000)

lats = []
lngs = []
for it in range(0, len(ts)):
    ta = func.get_ta(a, e, cons.mu, t0, ts[it])
    r = elements_to_vec(a, e, i, omega, w, ta)
    lng, lat = get_geo_coords(r, ts[it])

    lats.append(np.degrees(lat))
    lngs.append(np.degrees(lng))

ax = plt.axes(projection=ccrs.PlateCarree())
xticks = np.arange(-180, 180, 10)
yticks = np.arange(-90, 90, 10)
ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)

ax.set_global()
ax.coastlines()

plt.plot(lngs, lats, transform=ccrs.Geodetic())
plt.scatter(lngs[-1], lats[-1], marker='x', transform=ccrs.Geodetic(), color='r', s=50)
plt.show()
