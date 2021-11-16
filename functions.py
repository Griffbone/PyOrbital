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

    if e <= 1:
        p = a*(1 - e**2)
    else:
        # p = a*(e**2 - 1)
        p = a*(1 - e**2)
        asm = np.arccos(-1/e)
        thetas = np.linspace(-asm + abs(margin), asm - abs(margin), len(thetas))

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


def elements_to_eci_vel(a, e, i, lan, w, ta, mu=cns.mu):
    """ Function to calculate ECI state vector from orbital elements
            :param a: semimajor axis (m)
            :param e: eccentricity
            :param i: inclination (deg)
            :param lan: longnitude of ascending node (deg)
            :param w: argument of periapsis
            :param ta: true anomaly (deg)
            :param mu: central body gravitational parameter

            :return vx: x component of velocity vector
            :return vy: y component of velocity vector
            :return vz: z component of velocity vector
        """

    if e <= 1:
        h = np.sqrt(mu * a * (1 + e ** 2))
    else:
        h = np.sqrt(mu*a*(e**2 - 1))
        asm = np.degrees(np.arccos(-1 / e))

        if asm <= wrap_to_360(ta) <= 360 - asm:
            raise ValueError('Cannot find ECI velocity at specified true anomaly; orbit is hyperbolic.')

    vx = (mu / h) * (-np.sin(np.radians(ta)))
    vy = (mu / h) * (e + np.cos(np.radians(ta)))

    vx, vy, vz = perifocal_to_eci(lan, i, w, vx, vy)

    return vx, vy, vz


def elements_to_eci_pos(a, e, i, lan, w, ta):
    """ Function to calculate ECI state vector from orbital elements
            :param a: semimajor axis (m)
            :param e: eccentricity
            :param i: inclination (deg)
            :param lan: longnitude of ascending node (deg)
            :param w: argument of periapsis
            :param ta: true anomaly (deg)
            :param mu: central body gravitational parameter

            :return x: x component of position vector
            :return y: y component of position vector
            :return z: z component of position vector
        """

    if e <= 1:
        p = a*(1 - e**2)
    else:
        p = a*(e**2 - 1)
        asm = np.degrees(np.arccos(-1/e))

        if asm <= wrap_to_360(ta) <= 360 - asm:
            raise ValueError('Cannot find ECI coordinates at specified true anomaly; orbit is hyperbolic.')

    r = p / (1 + e * np.cos(np.radians(ta)))

    x = r * np.cos(np.radians(ta))
    y = r * np.sin(np.radians(ta))

    x, y, z = perifocal_to_eci(lan, i, w, x, y)

    return x, y, z


def vector_to_elements(r, v, mu):
    """ Function to determine orbital elements from state vector
        :param r: position vector
        :param v: velocity vector
        :param mu: central body gravitational parameter

        :return a: semimajor axis
        :return e: eccentricity
        :return i: inclination (deg)
        :return lan: longitude of ascending node (deg)
        :return w: argument of periapsis (deg)
        :return ta: true anomaly (deg)
    """

    rvec = r
    r = np.linalg.norm(rvec)

    vvec = v
    v = np.linalg.norm(vvec)

    vr = np.dot(rvec, vvec)/r

    # Angular momentum
    hvec = np.cross(rvec, vvec)
    h = np.linalg.norm(hvec)

    # Inclination
    i = np.arccos(hvec[2]/h)

    # LAN
    nvec = np.cross([0, 0, 1], hvec)
    n = np.linalg.norm(nvec)

    if n == 0:
        lan = 0
    else:
        lan = np.arccos(nvec[0] / n)
        if nvec[1] < 0:
            lan = 2*np.pi - lan

    # eccentricity
    evec = (1/mu)*((v**2 - mu/r)*rvec - r*vr*vvec)
    e = np.linalg.norm(evec)

    # argument of perigee
    if n == 0:
        # if orbit is in the equitorial plane, argument of periapsis is given from the equinox direction
        w = vector_angle(evec/e, np.array([1, 0, 0]))
    elif e > 0:
        w = np.arccos(np.dot(nvec, evec)/(n*e))
        if evec[2] < 0:
            w = 2*np.pi - w
    else:
        w = 0

    # True anomaly
    if e == 0:
        cp = np.cross(nvec, rvec)
        if cp[3] >= 0:
            ta = np.arccos(np.dot(nvec, rvec)/(n*r))
        else:
            ta = 2*np.pi - np.arccos(np.dot(nvec, rvec)/(n*r))
    else:
        ta = np.arccos(np.dot(evec, rvec)/(e*r))
        if vr < 0:
            ta = 2*np.pi - ta

    # semimajor axis
    a = 1/(2/r - v**2/mu)

    if e > 1:
        a = -a

    return a, e, np.degrees(i), np.degrees(lan), np.degrees(w), np.degrees(ta)


# def vector_to_elements(r, vel, mu):
#     """ Function to calculate orbital elmements given ECI state vector
#         :param r: position vector
#         :param vel: velocity vector
#         :param mu: gravitational parameter
#
#         :return a: semimajor axis
#         :return e: eccentricity
#         :return i: inclination
#         :return omega: longnitude of ascending node
#         :return w: argument of periapsis
#         :return M: mean anomaly
#     """
#
#     h = np.cross(r, vel)
#     evec = np.cross(vel, h)/mu - r/np.linalg.norm(r)
#     n = np.transpose(np.array([-h[1], h[0], 0]))
#
#     if np.dot(r, vel) > 0:
#         v = np.arccos(np.dot(evec, r)/(np.linalg.norm(evec)*np.linalg.norm(r)))
#     else:
#         v = 2*np.pi - np.arccos(np.dot(evec, r)/(np.linalg.norm(evec)*np.linalg.norm(r)))
#
#     # print(v)
#
#     i = np.arccos(h[2]/np.linalg.norm(h))
#     e = np.linalg.norm(evec)
#
#     if e < 1:
#         E = 2*np.arctan(np.tan(v/2)*((1 + e)/(1 - e))**(-1/2))
#         M = E - e * np.sin(E)
#     elif e > 1:
#         H = 2*np.arctanh(np.tan(v/2)*((1 + e)/(e - 1))**(-1/2))
#         M = e*np.sinh(H) - H
#     elif e == 1:
#         raise ValueError("Well, you're fucked")
#
#     if i == 0 or i == np.pi:
#         omega = 0
#     elif n[1] >= 0:
#         omega = np.arccos(n[0]/np.linalg.norm(n))
#     else:
#         omega = 2*np.pi - np.arccos(n[0]/np.linalg.norm(n))
#
#     if i == 0 or i == np.pi:
#         x = np.array([1, 0, 0])
#         w = vector_angle(x, evec)
#     elif evec[2] >= 0:
#         w = np.arccos(np.dot(n, evec)/(np.linalg.norm(n)*np.linalg.norm(evec)))
#     else:
#         w = 2*np.pi - np.arccos(np.dot(n, evec)/(np.linalg.norm(n)*np.linalg.norm(evec)))
#
#     # M = E - e*np.sin(E)
#     a = 1/((2/np.linalg.norm(r)) - (np.linalg.norm(vel)**2/mu))
#
#     return a, e, i, omega, w, kepler_ta(e, M)


# ====================== MISCELLANEOUS ORBIT MATH ======================

def period(a, mu):
    """ Function to get orbital period from semimajor axis
        :param a: semimajor axis (m)
        :param mu: parent gravitational parameter

        :return T: period (seconds)
    """

    T = np.sqrt((4*np.pi**2*a**3)/mu)

    return T