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

    _, _, x, y = func.perifocal_coords(a, e, tas)
    z = np.zeros(len(x))

    xps = []
    yps = []
    zps = []

    for x, y, z in zip(x, y, z):
        r = np.array([x, y, z])
        xp, yp, zp = func.perifocal_to_eci(omega, i, w, r)

        xps.append(xp)
        yps.append(yp)
        zps.append(zp)

    return xps, yps, zps


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
        t = t0 + ts[i]/(60**2 * 24)
        phi, lam, _ = func.eci_to_lla(xs[i], ys[i], zs[i], t)

        phis.append(phi)
        lams.append(lam)

    return phis, lams

