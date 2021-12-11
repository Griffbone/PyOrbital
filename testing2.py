import numpy as np
import functions as func
import constants as cns
import matplotlib.pyplot as plt


def circle(r, cx, cy, style):
    ts = np.linspace(0, 2*np.pi, 100)
    xs = r*np.cos(ts) + cx
    ys = r*np.sin(ts) + cy

    plt.plot(xs, ys, style)


def vector(cx, cy, vec, scale=1, style='ko-'):
    plt.plot([cx, cx + vec[0]*scale], [cy, cy + vec[1]*scale], style, ms=3, markevery=[-1])


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

    # Probably a good idea to add checks to determine if inputs are reasonable
    #   hyperbolic ta not within asymptotes
    #   LAN is zero if inclination is zero

    tol = 1e-6

    if not (1 - tol <= e <= 1 + tol):
        p = a*(1 - e**2)
    else:
        raise ValueError("Well, You're fucked. This function cannot handle parabolic orbits")

    ta = np.radians(ta)

    rpf = np.array([p*np.cos(ta)/(1 + e*np.cos(ta)), p*np.sin(ta)/(1 + e*np.cos(ta)), 0])
    vpf = np.array([-np.sqrt(mu/p)*np.sin(ta), np.sqrt(mu/p)*(e + np.cos(ta)), 0])

    rotmat = rotz(-w) @ rotx(-i) @ rotz(-lan) #rotz(-lan) @ rotx(-i) @ rotz(-w)

    r = rotmat @ rpf
    v = rotmat @ vpf

    return r, v


r = np.array([1, 0, 1])

print((r - 1e-6 <= r).all() and (r <= r + 1e-6).all())
