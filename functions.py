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


def rotx(t):
    t = np.radians(t)
    cost = np.cos(t)
    sint = np.sin(t)

    return np.array([[1, 0, 0],
                     [0, cost, -sint],
                     [0, sint, cost]])


def rotz(t):
    t = np.radians(t)
    cost = np.cos(t)
    sint = np.sin(t)

    return np.array([[cost, -sint, 0],
                     [sint, cost, 0],
                     [0, 0, 1]])

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

def elements_to_vector(a, e, i, lan, w, ta, truelon, arglat, lonper, mu):
    """ Function to convert orbital elements to position and velocity
        :param a: semimajor axis
        :param e: eccentricity
        :param i: inclination
        :param lan: longitude of ascending node
        :param w: argument of periapsis
        :param ta: true anomaly
        :param truelon: true longitude
        :param arglat: argument of latitude
        :param lonper: longitude of periapsis
        :param mu: gravitational parameter

        :return r: position vector
        :return v: velocity vector
    """

    # Probably a good idea to add checks to determine if inputs are reasonable
    #   hyperbolic ta not within asymptotes
    #   LAN is zero if inclination is zero

    tol = 1e-6
    w = np.radians(w)
    lan = np.radians(lan)
    i = np.radians(i)

    if not (1 - tol <= e <= 1 + tol):
        p = a*(1 - e**2)
    else:
        raise ValueError("Well, You're fucked. This function cannot handle parabolic orbits")

    ta = np.radians(ta)

    rpf = np.array([p*np.cos(ta)/(1 + e*np.cos(ta)), p*np.sin(ta)/(1 + e*np.cos(ta)), 0])
    vpf = np.array([-np.sqrt(mu/p)*np.sin(ta), np.sqrt(mu/p)*(e + np.cos(ta)), 0])

    clan = np.cos(lan)
    slan = np.sin(lan)
    cosw = np.cos(w)
    sinw = np.sin(w)
    cosi = np.cos(i)
    sini = np.sin(i)

    rotmat = np.array([[clan*cosw - slan*sinw*cosi, -clan*sinw - slan*cosw*cosi, slan*sini],
                       [slan*cosw + clan*sinw*cosi, -slan*sinw + clan*cosw*cosi, -clan*sini],
                       [sinw*sini, cosw*sini, cosi]])

    r = rotmat @ rpf
    v = rotmat @ vpf

    return r, v


def vector_to_elements(rvec, vvec, mu):
    """ Function to convert position and velocity to classical orbital elements
        :param rvec: position vector
        :param vvec: velocity vector
        :param mu: gravitational parameter

        :return a: semimajor axis
        :return e: eccentricity
        :return i: inclination (deg)
        :return lan: longitude of ascending node (deg)
        :return w: argument of periapsis (deg)
        :return ta: true anomaly (deg)
        :return truelon: true longitude (deg)
        :return arglat: argument of latitude (deg)
        :return lonper: longitude of periapsis (deg)
    """

    tol = 1e-6
    r = np.linalg.norm(rvec)
    v = np.linalg.norm(vvec)
    E = v**2/2 - mu/r

    hvec = np.cross(rvec, vvec)
    h = np.linalg.norm(hvec)

    # Radial orbit
    if h < tol:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    nvec = np.cross([0, 0, 1], hvec)
    n = np.linalg.norm(nvec)

    evec = (1/mu)*((v**2 - mu/r)*rvec - np.dot(rvec, vvec)*vvec)
    e = np.linalg.norm(evec)

    cosi = hvec[2]/h
    i = np.degrees(np.arccos(cosi))

    # Semimajor axis
    if abs(E) > tol:
        # Elliptical/hyperbolic orbit
        a = -mu/(2*E)
    else:
        # Parabolic orbit
        a = np.inf

    # Semiminor axis
    if n > tol:
        # Inclined orbit case
        lan = np.degrees(np.arccos(nvec[0]/n))

        if nvec[1] < 0:
            lan = 360 - lan
    else:
        # Equatorial orbit case
        lan = np.nan

    # Argument of periapsis
    if (e > tol) and (i > tol):
        # Eccentric inclined case
        w = np.degrees(np.arccos(np.dot(nvec, evec)/(n*e)))

        if evec[2] < 0:
            w = 360 - w
    else:
        # Circular or non-inclined case
        w = np.nan

    # True anomaly
    if e > tol:
        # Eccentric orbit case
        ta = np.degrees(np.arccos(np.dot(evec, rvec)/(e*r)))

        if np.dot(rvec, vvec) < 0:
            ta = 360 - ta
    else:
        # Circular orbit case
        ta = np.nan

    # Argument of latitude (inclined orbit)
    if n > tol:   # (e <= tol) and (n > tol):
        arglat = np.degrees(np.arccos(np.dot(nvec, rvec)/(n*r)))

        if rvec[2] < 0:
            arglat = 360 - arglat
    else:
        arglat = np.nan

    # True longitude (circular equatorial orbit)
    if (e <= tol) and (n <= tol):
        truelon = np.degrees(np.arccos(rvec[0]/r))

        if (rvec[1] < 0) or (i > 90):
            truelon = 360 - truelon
    else:
        truelon = np.nan

    # True longitude of periapsis (elliptical equatorial orbit)
    if (e > tol) and (n < tol):
        lonper = np.degrees(np.arccos(evec[0] / e))

        if evec[1] < 0:
            lonper = 360 - lonper
    else:
        lonper = np.nan

    return a, e, i, lan, w, ta, arglat, truelon, lonper


# ====================== MISCELLANEOUS ORBIT MATH ======================

def period(a, mu):
    """ Function to get orbital period from semimajor axis
        :param a: semimajor axis (m)
        :param mu: parent gravitational parameter

        :return T: period (seconds)
    """

    T = np.sqrt((4*np.pi**2*a**3)/mu)

    return T