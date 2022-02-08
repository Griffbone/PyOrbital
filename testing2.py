import numpy as np
import constants as cns
import astrotime as astt

a = cns.smae
f = 1/cns.flat
b = a*(1 - f)


404-385-3314
404-650-5689


def lla_to_ecef(lat, lon, alt):
    """ Function to convert geodetic coordinates and altitude to ECEF coordinates
        :param lat: geodetic latitude
        :param lon: geodetic longitude
        :param alt: altitude above ellipsoid (m)
    """
    e = 1/cns.flat
    lat = np.radians(lat)
    lon = np.radians(lon)

    xp = (cns.smae / np.sqrt(1 - e**2 * np.sin(lat)**2) + alt)*np.cos(lat)
    z = (cns.smae * (1 - e**2) / np.sqrt(1 - e**2 * np.sin(lat)**2) + alt)*np.sin(lat)

    x = xp*np.cos(lon)
    y = xp*np.sin(lon)

    return x, y, z


def lla_to_eci(lat, lon, alt, jdn):
    """ Function to convert geodetic coordinates and altitude to ECEF coordinates
        :param lat: geodetic latitude
        :param lon: geodetic longitude
        :param alt: altitude above ellipsoid (m)
        :param jdn: Julian day number
    """
    e = 1/cns.flat
    lat = np.radians(lat)
    lon = astt.theta_g(jdn) + np.radians(lon)

    xp = (cns.smae / np.sqrt(1 - e**2 * np.sin(lat)**2) + alt)*np.cos(lat)
    z = (cns.smae * (1 - e**2) / np.sqrt(1 - e**2 * np.sin(lat)**2) + alt)*np.sin(lat)

    x = xp*np.cos(lon)
    y = xp*np.sin(lon)

    return x, y, z


jdn = astt.date_to_jd(1, 0, 2016, 6, 36, 25)
print(np.radians(astt.theta_g_2(jdn)))


print(lla_to_eci(0, -57.296, 6.378e3, astt.date_to_jd(1, 2, 1970, 6, 0, 0)))
