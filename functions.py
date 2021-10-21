import numpy as np
from scipy.optimize import fsolve
import astrotime as at


def wrap_to_360(theta):
    """ Function to wrap an angle to [0, 360]
        :param theta: some angle (deg)
        :return theta: wrapped angle
    """

    while theta > 360:
        theta -= 360

    while theta < 0:
        theta += 360

    return theta


def rotmat(t, axis):
    """Function to create a rotation matrix about an axis
        :param t: angular displacement in RADIANS
        :type t: float
        :param axis: axis of rotation
        :type axis: str
    """

    if axis.upper() == 'X':
        return np.array([
            [1, 0, 0],
            [0, np.cos(t), -np.sin(t)],
            [0, np.sin(t), np.cos(t)]
        ])

    elif axis.upper() == 'Y':
        return np.array([
            [np.cos(t), 0, np.sin(t)],
            [0, 1, 0],
            [-np.sin(t), 0, np.cos(t)]
        ])

    elif axis.upper() == 'Z':
        return np.array([
            [np.cos(t), -np.sin(t), 0],
            [np.sin(t), np.cos(t), 0],
            [0, 0, 1]
        ])


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


def perifocal_to_eci(omega, i, w, x, y):
    """" Function to transform perifocal position vector to ECI position vector
        :param omega: right ascension of ascending node (deg)
        :param i: inclination (deg)
        :param w: argument of periapsis (deg)
        :param r: 1x3 perifocal position vector

        :return r_eci: position vector in the ECI frame
    """

    Z1 = rotmat(np.radians(omega), 'z')
    X1 = rotmat(np.radians(i), 'x')
    Z2 = rotmat(np.radians(w), 'z')

    r = np.array([x, y, 0])

    R = Z1 @ X1 @ Z2
    r_eci = R @ r

    return r_eci[0], r_eci[1], r_eci[2]


def elements_to_vector(a, e, i, lan, w, ta):
    """ Function to calculate ECI state vector from orbital elements
        :param a: semimajor axis (m)
        :param e: eccentricity
        :param i: inclination (deg)
        :param lan: longnitude of ascending node (deg)
        :param w: argument of periapsis
        :param ta: true anomaly (deg)

        :return x: x component of state vector
        :return y: y component of state vector
        :return z: z component of state vector
    """

    p = (a * np.sqrt(1 - e ** 2)) ** 2 / a
    r = p / (1 + e * np.cos(np.radians(ta)))

    x = r*np.cos(np.radians(ta))
    y = r*np.sin(np.radians(ta))

    x, y, z = perifocal_to_eci(lan, i, w, x, y)

    return x, y, z


def kepler_ta(e, ma):
    """ Solve Kepler's equation for Eccentric anomaly and True Anomaly
        :param e: eccentricity
        :param ma: mean anomaly (deg)

        :return TA: true anomaly (deg)
    """

    ma = np.radians(ma)

    E = fsolve(lambda E: E - e*np.sin(E) - ma, np.array([0]))

    if len(E) > 1:
        E = E[0]

    TA = np.degrees(2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2)))
    TA = wrap_to_360(TA)

    return TA[0]


def eci_to_ll(x, y, z, jdn):
    """ Function to convert ECI coordinates to lat/lon
        :param x: ECI x position (m)
        :param y: ECI y position (m)
        :param z: ECI z position (m)
        :param jdn: julian day number

        :return phi: latitude (deg)
        :return lam: longitude (deg)
    """

    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arcsin(z/r)

    lam = (np.degrees(np.arctan2(y, x)) - at.theta_g(jdn))
    lam = wrap_to_360(lam)

    return np.degrees(phi), lam


def ta_to_ma(e, ta):
    """ Functtion to get mean anomaly from true anomaly
        :param ta: true anomaly (deg)
        :return ma: mean anomaly (deg)
    """
    ta = np.radians(ta)
    ea = np.arctan(np.tan(ta/2)/np.sqrt((1+e)/(1-e)))*2
    ma = wrap_to_360(np.degrees(ea - e*np.sin(ea)))

    return ma


def vector_to_elements(r, vel, mu):
    h = np.cross(r, vel)
    evec = np.cross(vel, h)/mu - r/np.linalg.norm(r)
    # n = np.cross(np.transpose([0, 0, 1]), h)
    n = np.transpose(np.array([-h[1], h[0], 0]))

    if np.dot(r, vel) > 0:
        v = np.arccos(np.dot(evec, r)/(np.linalg.norm(evec)*np.linalg.norm(r)))
    else:
        v = 2*np.pi - np.arccos(np.dot(evec, r)/(np.linalg.norm(evec)*np.linalg.norm(r)))

    i = np.arccos(h[2]/np.linalg.norm(h))
    e = np.linalg.norm(evec)
    E = 2*np.arctan(np.tan(v/2)/np.sqrt((1+2)/(1-e)))

    if n[1] >= 0:
        omega = np.arccos(n[0]/np.linalg.norm(n))
    else:
        omega = 2*np.pi - np.arccos(n[0]/np.linalg.norm(n))

    if evec[2] >= 0:
        w = np.arccos(np.dot(n, evec)/(np.linalg.norm(n)*np.linalg.norm(evec)))
    else:
        w = 2*np.pi - np.arccos(np.dot(n, evec)/(np.linalg.norm(n)*np.linalg.norm(evec)))

    M = E - e*np.sin(E)
    a = 1/((2/np.linalg.norm(r)) - (np.linalg.norm(vel)**2/mu))

    return a, e, i, omega, w, M
