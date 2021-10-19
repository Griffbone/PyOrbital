import numpy as np
from scipy.optimize import fsolve


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


def kepler_ta(e, ma):
    """ Solve Kepler's equation for Eccentric anomaly and True Anomaly
        :param e: eccentricity
        :param ma: mean anomaly
        :return TA: true anomaly
    """

    E = fsolve(lambda E: E - e*np.sin(E) - ma, np.array([0]))

    if len(E) > 1:
        E = E[0]

    TA = np.arctan(np.sqrt((1 + e)/(1 - e))*np.tan(E/2))*2

    while TA < 0:
        TA += np.pi*2

    while TA > 2*np.pi:
        TA -= 2*np.pi

    return TA[0]


def rot_orbit(omega, i, w, xs, ys, zs):
    Z1 = rotmat(omega, 'z')
    X1 = rotmat(i, 'x')
    Z2 = rotmat(w, 'z')

    R = Z1 @ X1 @ Z2

    # rotate orbit in IJK space using rotation matrix
    xps = []
    yps = []
    zps = []

    for x,y,z in zip(xs, ys, zs):
        vnew = R@np.array([x, y, z])

        xps.append(vnew[0])
        yps.append(vnew[1])
        zps.append(vnew[2])

    return xps, yps, zps


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
