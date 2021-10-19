import numpy as np
import functions as func


def perifocal_coords(a, e, thetas=np.linspace(0, 2*np.pi, 100)):
    """ Function to return perifocal polar coordinates of an orbit
        :param a: semimajor axis
        :param e: eccentricity
        :param thetas: array of true anomalies (rad)

        :return rs: array of distances
        :return thetas: array of true anomalies
        :return xs: array of x locations
        :return ys: array of y locations
    """

    p = (a*np.sqrt(1 - e**2))**2 / a
    rs = p/(1 + e*np.cos(thetas))

    xs = rs*np.cos(thetas)
    ys = rs*np.sin(thetas)

    return rs, thetas, xs, ys


def perifocal_to_eci(omega, i, w, r):
    """" Function to transform perifocal position vector to ECI position vector
        :param omega: right ascension of ascending node (deg)
        :param i: inclination (deg)
        :param w: argument of periapsis (deg)
        :param r: 1x3 perifocal position vector

        :return r_eci: position vector in the ECI frame
    """

    Z1 = func.rotmat(np.radians(omega), 'z')
    X1 = func.rotmat(np.radians(i), 'x')
    Z2 = func.rotmat(np.radians(w), 'z')

    R = Z1 @ X1 @ Z2
    r_eci = R @ r

    return r_eci


def plot_orbit(a, e, i, omega, w, tas=np.linspace(0, 2*np.pi, 100)):
    """ Function to plot an orbit from the orbit elements
        :param a: semimajor axis
        :param e: eccentricity
        :param i: inclination
        :param w: argument of periapsis
        :param tas: true anomalies to plot at

        :return x: vector of x positions
        :return y: vector of y positions
        :return z: vector of z positions
    """

    _, _, x, y = perifocal_coords(a, e, tas)
    z = np.zeros(len(x))

    xps = []
    yps = []
    zps = []

    for x, y, z in zip(x, y, z):
        r = np.array([x, y, z])
        xp, yp, zp = perifocal_to_eci(omega, i, w, r)

        xps.append(xp)
        yps.append(yp)
        zps.append(zp)

    return xps, yps, zps

