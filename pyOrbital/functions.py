import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.optimize import fsolve


def get_ta(a, e, mu, t0, t):
    """ Function to get true anomaly from orbital elements and time
        :param a: semimajor axis (m)
        :param e: eccentricity
        :param mu: gravitational parameter
        :param t0: time of periapsis passage (JD)
        :param t: current time (JD)
    """

    T = 2*np.pi*np.sqrt(a**3/mu)
    t = (t - t0)*24*60*60

    n = 2*np.pi/T
    M = n*t
    E = fsolve(lambda E: E - e*np.sin(E) - M, np.array([M]))

    if len(E) > 1:
        E = E[0]

    TA = np.arctan(np.sqrt((1 + e)/(1 - e))*np.tan(E/2))*2

    while TA < 0:
        TA += np.pi*2

    while TA > 2*np.pi:
        TA -= 2*np.pi

    return TA[0]


def plot_orbit_in_plane(a, e, t0, t1, t2, mu):
    b = a*np.sqrt(1 - e**2)
    p = b**2/a

    # ta1 = get_ta(a, e, mu, t0, t1)
    # ta2 = get_ta(a, e, mu, t0, t2)
    ta1 = 0
    ta2 = np.pi*2

    thetas = np.linspace(ta1, ta2, 100)
    r = p/(1 - e*np.cos(thetas))

    x = -r*np.cos(thetas)
    y = r*np.sin(thetas)

    return x, y, np.zeros(len(x))


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


def plot_planet():
    ts = np.linspace(0, 2*np.pi, 1000)
    rs = np.ones(1000)*6371
    plt.polar(ts, rs)


# def elements_to_vector(a, e, i, omega, w, E):


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


class orbit_from_elements():
    def __init__(self, a, e, i, lan, argPe, ma=0, ta=0):
        """Generate an orbit object from orbital elements
            :param a: semimajor axis
            :type a: float
            :param e: orbit eccentricity
            :type e: float
            :param i: orbit inclination
            :type i: float
            :param lan: longnitude of the ascending node
            :type lan: float
            :param argPe: argument of periapsis
            :type argPe: float
            :param ma: mean anomaly
            :type ma: float
            :param ta: true anomaly
            :type ta: float
        """

        self.e = e
        self.i = i
        self.lan = lan
        self.argPe = argPe
        self.ma = ma
        self.ta = ta

        self.b = a*np.sqrt(1 - e**2)    # semiminor axis
        self.p = self.b**2/a            # parameter

    def plot(self, plottype='inplane', style='-', lineweight=2, n=100):
        """Function to plot an orbit using matplotlib
            :param plottype: plot type (inplane/)
            :type plottype: str
            :param style: plot style
            :type style: str
        """

        if type(plottype).__name__ != 'str':
            print('ERROR: plottype type error; expected str but got ' + type(plottype).__name__)
            return None

        thetas = np.linspace(0, np.pi*2, n)
        rs = self.p/(1 - self.e*np.cos(thetas))

        xs = -rs*np.cos(thetas)
        ys = rs*np.sin(thetas)
        zs = np.zeros(n)

        xs, ys, zs = rot_orbit(self.lan, self.i, self.argPe, xs, ys, zs)

        return xs, ys, zs

        # if plottype.upper() == 'INPLANE':
        #     plt.plot(np.array(xs)/1000, np.array(ys)/1000, style)
        #     plt.gca().set_aspect('equal')
        #     plt.gca().add_patch(patches.Circle((0, 0), 6371, alpha=0.5, edgecolor='k'))
        #     plt.xlabel('X [km]')
        #     plt.ylabel('Y [km]')
        #     plt.title('Cartesian Orbit Plot (In-Plane)')
        #     plt.show()
        # elif plottype.upper() == 'PLANES':
        #     plt.subplot(1,3,1)
        #     plt.plot(np.array(xs)/1000, np.array(ys)/1000, style)
        #     plt.gca().set_aspect('equal')
        #     plt.gca().add_patch(patches.Circle((0, 0), 6371, alpha=0.5, edgecolor='k'))
        #
        #     plt.subplot(1,3,2)
        #     plt.plot(np.array(xs)/1000, np.array(zs)/1000, style)
        #     plt.gca().set_aspect('equal')
        #     plt.gca().add_patch(patches.Circle((0, 0), 6371, alpha=0.5, edgecolor='k'))
        #
        #     plt.subplot(1,3,3)
        #     plt.plot(np.array(ys)/1000, np.array(zs)/1000, style)
        #     plt.gca().set_aspect('equal')
        #     plt.gca().add_patch(patches.Circle((0, 0), 6371, alpha=0.5, edgecolor='k'))
        #
        #     plt.show()
        # else:
        #     print('ERROR: unrecognized plot type: ' + plottype)
        #     return None
