import math

import prediction as pred
import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import plotting as aplt
import functions as func


def stumpff_c_prime(z):
    c = 0
    n = 0

    while n <= 50:
        c += ((-1)**(n+1))*(n+1)*(z**n)/math.factorial(2*n + 4)
        n += 1

    return c


def stumpff_s_prime(z):
    s = 0
    n = 0

    while n <= 50:
        s += ((-1)**(n+1))*(n+1)*(z**n)/math.factorial(2*n + 5)
        n += 1

    return s


def lambert_pit(r1, r2, tof, mu):
    """ Function to solve Lambert's problem via P-iteration (BMW)
        :param r1: initial position vector
        :param r2: final position vector
        :param tof: time of flight
        :param mu: gravitational parameter

        :return v1: initial velocity
        :return v2: final velocity
    """
    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)



def lambert_uv(r1, r2, tof, mu, dir=1):
    """ Function to solve Gauss' problem via universal variables
        :param r1: initial position vector
        :param r2: final position vector
        :param tof: time of flight
        :param mu: gravitational parameter
        :param dir: long way or short way (1 = long, 2 = short)

        :return v1: initial velocity
        :return v2: final velocity
    """

    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)
    d_theta = np.arccos(np.dot(r1, r2)/(r1_norm*r2_norm))

    if dir == 1:
        d_theta = 2*np.pi - d_theta
    elif dir == 2:
        d_theta = d_theta

    # Handle error of being on a degenerate conic
    if d_theta == 0:
        return None

    A = (np.sqrt(r1_norm*r2_norm)*np.sin(d_theta))/(np.sqrt(1 - np.cos(d_theta)))

    # Determine Z
    z = 0
    t = 2*tof
    n = 0

    while abs(t - tof) > 1e-4:
        c = pred.stumpff_c(z)
        s = pred.stumpff_s(z)
        y = r1_norm + r2_norm - A*(1 - z*s)/np.sqrt(c)

        # Handle error of y being less than one
        if y < 1:
            return None

        x = np.sqrt(y/c)
        t = (1/np.sqrt(mu))*(x**3*s + A*np.sqrt(y))

        sprime = stumpff_s_prime(z)
        cprime = stumpff_c_prime(z)
        dtdz = (1/np.sqrt(mu))*((x**3)*(sprime - (3*s*cprime)/(2*c)) + (A/8)*((3*s*np.sqrt(y))/c + A/x))
        z = z + (tof - t)/dtdz

        n += 1

    c = pred.stumpff_c(z)
    s = pred.stumpff_s(z)
    y = r1_norm + r2_norm - A * (1 - z * s) / np.sqrt(c)

    f = 1 - y/r1_norm
    g = A*np.sqrt(y/mu)
    gdot = 1 - y/r2_norm

    v1 = (r2 - f*r1)/g
    v2 = (gdot*r2 - r1)/g

    return v1, v2
