import numpy as np
import functions as func


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

