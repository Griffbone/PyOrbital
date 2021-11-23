import numpy as np
from scipy.optimize import fsolve
import astrotime as at
import constants as cns
import matplotlib.pyplot as plt

# ====================== ELEMENT SET CLASS ======================


class Elements:
    def __init__(self, a, e, i, lan, w, ta, mu):
        self.a = a
        self.e = e
        self.i = i
        self.lan = lan
        self.w = w
        self.ta = ta
        self.mu = mu

        self.ap = a*(1 + e)
        self.pe = a*(1 - e)

        self.T = period(abs(a), mu)
        self.n = 2*np.pi/self.T

    def perifocal_plot(self, n=1000):
        _, _, x, y = perifocal_coords(self.a, self.e, np.linspace(0, 2*np.pi, n))
        plt.plot(x, y)
        plt.axis('equal')

    def get_eci(self, ta=None):
        if ta is None:
            ta = self.ta

        x, y, z = elements_to_eci_pos(self.a, self.e, self.i, self.lan, self.w, ta)

        return x, y, z


# ====================== BASIC MATH FUNCTIONS ======================

def vector_angle(v1, v2):
    """ Fuction to calculate an angle between two vectors
        :param v1: first vector
        :param v2: second vector

        :return theta: angle (rad)
    """

    dot = np.dot(v1, v2)
    mag = np.linalg.norm(v1)*np.linalg.norm(v2)

    if dot == 0:
        theta = 0
    else:
        theta = np.arccos(dot/mag)

    return theta


def wrap_to_360(theta):
    """ Function to wrap an angle to [0, 360]
        :param theta: some angle (deg)
        :return theta: wrapped angle (deg)
    """

    while theta > 360:
        theta -= 360

    while theta < 0:
        theta += 360

    return theta


def wrap_to_2pi(theta):
    """ Function to wrap an angle to [0, 2*pi]
        :param theta: some angle (rad)
        :return theta: wrapped angle (rad)
    """

    while theta > 2*np.pi:
        theta -= 2*np.pi

    while theta < 0:
        theta += 2*np.pi

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

# ====================== ORBITAL POSITION FUNCTIONS ======================


def perifocal_coords(a, e, thetas=np.linspace(0, 2*np.pi, 100), margin=0.1):
    """ Function to return perifocal polar coordinates of an orbit
        :param a: semimajor axis
        :param e: eccentricity
        :param thetas: array of true anomalies (rad)
        :param margin: margin to add to asymptotes for a hyperbolic orbit

        :return rs: array of distances
        :return thetas: array of true anomalies
        :return xs: array of x locations
        :return ys: array of y locations
    """

    p = a*(1 - e**2)
    #
    # if e > 1:
    #     p = a*(1 - e**2)

    # add some logic for hyperbolic orbits asymptotes

        # asm = np.arccos(-1/e)
        # thetas = np.linspace(-asm + abs(margin), asm - abs(margin), len(thetas))

    rs = p/(1 + e*np.cos(thetas))

    xs = rs*np.cos(thetas)
    ys = rs*np.sin(thetas)

    return rs, thetas, xs, ys


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


def ta_to_ma(e, ta):
    """ Functtion to get mean anomaly from true anomaly
        :param ta: true anomaly (deg)
        :return ma: mean anomaly (deg)
    """
    ta = np.radians(ta)
    ea = np.arctan(np.tan(ta/2)/np.sqrt((1+e)/(1-e)))*2
    ma = wrap_to_360(np.degrees(ea - e*np.sin(ea)))

    return ma


def t_between(a, e, ta1, ta2, mu):
    """ Function to get time between two true anomalies
        :param ta1: first true anomaly (deg)
        :param ta2: second true anomaly (deg)
        :return Dt: time of flight
    """

    T = period(a, mu)
    n = 360/T

    m1 = ta_to_ma(e, ta1)
    m2 = ta_to_ma(e, ta2)

    t0 = m1/n
    t1 = m2/n

    Dt = t1 - t0

    if Dt < 0:
        Dt += T

    return Dt


def ta_change(a, e, ta1, Dt, mu):
    """ Function to get change in true anomaly from a change in time
        :param a: semimajor axis
        :param e: eccentricity
        :param ta1: original true anomaly (deg)
        :param Dt: change in time (s)
        :param mu: central gravitational parameter

        :return ta2: second true anomaly
    """
    T = period(a, mu)
    n = 360/T

    m1 = ta_to_ma(e, ta1)
    m2 = m1 + n*Dt

    ta2 = kepler_ta(e, m2)

    return ta2

# ====================== ORBITAL REFERENCE FRAME FUNCTIONS ======================


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


def eci_to_lla(x, y, z, jdn):
    """ Function to convert ECI coordinates to lat/lon
        ASSUMES SPHERICAL EARTH (for now)

        :param x: ECI x position (m)
        :param y: ECI y position (m)
        :param z: ECI z position (m)
        :param jdn: julian day number

        :return phi: latitude (deg)
        :return lam: longitude (deg)
        :return a: altitude (m)
    """

    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arcsin(z/r)

    lam = (np.degrees(np.arctan2(y, x)) - at.theta_g(jdn))
    lam = wrap_to_360(lam)

    a = r - cns.re

    return np.degrees(phi), lam, a


# ====================== ORBITAL ELEMENT CONVERSION ======================

def elements_to_vector(a, e, i, lan, w, ta, mu):
    """ Function to convert orbital elements to position and velocity
        :param a: semimajor axis
        :param e: eccentricity
        :param i: inclination
        :param lan: longitude of ascending node
        :param w: argument of periapsis
        :param ta: true anomaly
        :param mu: gravitational parameter

        :return r: position vector
        :return v: velocity vector
    """
    pass


def vector_to_elements(rvec, vvec, mu):
    """ Function to convert position and velocity to classical orbital elements

        DOUBLE CHECK TRUE ANOMALY SPECIAL CASES

        :param r: position vector
        :param v: velocity vector
        :param mu: gravitational parameter

        :return a: semimajor axis
        :return e: eccentricity
        :return i: inclination (deg)
        :return lan: longitude of ascending node (deg)
        :return w: argument of periapsis (deg)
        :return ta: true anomaly (deg)
    """
    tol = 1e-6

    r = np.linalg.norm(rvec)
    v = np.linalg.norm(vvec)
    E = v**2/2 - mu/r

    hvec = np.cross(rvec, vvec)
    h = np.linalg.norm(hvec)

    nvec = np.cross([0, 0, 1], hvec)
    n = np.linalg.norm(nvec)

    evec = (1/mu)*((v**2 - mu/r)*rvec - np.dot(rvec, vvec)*vvec)
    e = np.linalg.norm(evec)

    # Semimajor axis
    if (e > 1 - tol) and (e < 1 + tol):
        # Parabolic orbit case
        a = np.inf
    else:
        # Elliptic/hyperbolic case
        a = -mu/(2*E)

    # Inclination
    cosi = hvec[2]/h
    i = np.degrees(np.arccos(cosi))

    # LAN
    if i == 0 or i == 180:
        # equatorial orbit case
        lan = 0
    else:
        # general case
        coslan = nvec[0]/n
        lan = np.degrees(np.arccos(coslan))

        if nvec[1] < 0:
            lan = 360 - lan

    # Argument of periapsis
    if e <= tol:
        # circular orbit
        w = 0
    elif i == 0 or i == 180:
        # equatorial orbit
        cosw = evec[0]/e
        w = np.degrees(np.arccos(cosw))

        if evec[1] < 0:
            w = 360 - w
    else:
        # general case
        cosw = np.dot(nvec, evec)/(n*e)
        w = np.degrees(np.arccos(cosw))

        if evec[2] < 0:
            w = 360 - w

    # True anomaly
    if e <= tol and i != 0 and i != 180:
        # circular inclined orbit
        costa = np.dot(nvec, rvec)/(n*r)
        ta = np.degrees(np.arccos(costa))

        if rvec[2] < 0:
            ta = 360 - ta
    elif e <= tol and (i == 0 or i == 180):
        # circular equatorial orbit
        costa = rvec[0] / r
        ta = np.degrees(np.arccos(costa))

        if rvec[2] < 0:
            ta = 360 - ta
    else:
        # general case
        costa = np.dot(evec, rvec)/(e*r)
        ta = np.degrees(np.arccos(costa))

        if np.dot(rvec, vvec) < 0:
            ta = 360 - ta

    return a, e, i, lan, w, ta


# ====================== MISCELLANEOUS ORBIT MATH ======================

def period(a, mu):
    """ Function to get orbital period from semimajor axis
        :param a: semimajor axis (m)
        :param mu: parent gravitational parameter

        :return T: period (seconds)
    """

    T = np.sqrt((4*np.pi**2*a**3)/mu)

    return T