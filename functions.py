import numpy as np
from scipy.optimize import fsolve
import astrotime as at
import constants as cns
import matplotlib.pyplot as plt


# ====================== ELEMENT SET CLASS ======================
class Elements:
    pass


# ====================== BASIC MATH FUNCTIONS ======================
def wrap_to_360(theta):
    """ Function to wrap an angle to [0, 360]
        :param theta: some angle (deg)
        :return theta_w: wrapped angle (deg)
    """

    theta_w = theta % 360

    return theta_w


def wrap_to_2pi(theta):
    """ Function to wrap an angle to [0, 2*pi]
        :param theta: some angle (rad)
        :return theta_w: wrapped angle (rad)
    """

    theta_w = theta % (2*np.pi)

    return theta_w


def rot_1(theta):
    """ Function to create a counterclockwise rotation matrix about the X-axis
        :param theta: rotation angle (rad)
        :return rotmat: clockwise rotation matrix
    """

    ctheta = np.cos(theta)
    stheta = np.sin(theta)

    rotmat = np.array([[1, 0, 0], [0, ctheta, stheta], [0, -stheta, ctheta]])

    return rotmat


def rot_2(theta):
    """ Function to create a counterclockwise rotation matrix about the Y-axis
        :param theta: rotation angle (rad)
        :return rotmat: clockwise rotation matrix
    """

    ctheta = np.cos(theta)
    stheta = np.sin(theta)

    rotmat = np.array([[ctheta, 0, -stheta], [0, 1, 0], [stheta, 0, ctheta]])

    return rotmat


def rot_3(theta):
    """ Function to create a counterclockwise rotation matrix about the Z-axis
        :param theta: rotation angle (rad)
        :return rotmat: clockwise rotation matrix
    """

    ctheta = np.cos(theta)
    stheta = np.sin(theta)

    rotmat = np.array([[ctheta, stheta, 0], [-stheta, ctheta, 0], [0, 0, 1]])

    return rotmat


def vec_rot_1(theta, v):
    """ Perform a vector rotation about the X-axis.
        Functionally the same as using a rotation matrix output from rot_1, but has slightly better performance.

        :param theta: rotation angle (rad)
        :param v: vector to rotate
        :return vp: rotated vector
    """

    x = v[0]
    y = v[1]
    z = v[2]

    ctheta = np.cos(theta)
    stheta = np.sin(theta)

    yp = ctheta*y + stheta*z
    zp = -stheta*y + ctheta*z
    vp = np.array([x, yp, zp])

    return vp


def vec_rot_2(theta, v):
    """ Perform a vector rotation about the Y-axis.
        Functionally the same as using a rotation matrix output from rot_2, but has slightly better performance.

        :param theta: rotation angle (rad)
        :param v: vector to rotate
        :return vp: rotated vector
    """

    x = v[0]
    y = v[1]
    z = v[2]

    ctheta = np.cos(theta)
    stheta = np.sin(theta)

    xp = ctheta*x - stheta*z
    zp = stheta*x + ctheta*z
    vp = np.array([xp, y, zp])

    return vp


def vec_rot_3(theta, v):
    """ Perform a vector rotation about the Z-axis.
        Functionally the same as using a rotation matrix output from rot_3, but has slightly better performance.

        :param theta: rotation angle (rad)
        :param v: vector to rotate
        :return vp: rotated vector
    """

    x = v[0]
    y = v[1]
    z = v[2]

    ctheta = np.cos(theta)
    stheta = np.sin(theta)

    xp = ctheta*x + stheta*y
    yp = -stheta*x + ctheta*z
    vp = np.array([xp, yp, z])

    return vp


# ====================== ORBITAL POSITION FUNCTIONS ======================

# ====================== ORBITAL REFERENCE FRAME FUNCTIONS ======================
def sez_to_ecef():
    pass


def sez_to_eci():
    pass


def perifocal_to_eci():
    pass


def eci_to_lla():
    pass


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
