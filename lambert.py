import math

import prediction as pred
import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import plotting as aplt
import functions as func
import prediction as pred
import propagators


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
            return None, None

        x = np.sqrt(y/c)
        t = (1/np.sqrt(mu))*(x**3*s + A*np.sqrt(y))

        sprime = stumpff_s_prime(z)
        cprime = stumpff_c_prime(z)
        dtdz = (1/np.sqrt(mu))*((x**3)*(sprime - (3*s*cprime)/(2*c)) + (A/8)*((3*s*np.sqrt(y))/c + A/x))
        z = z + (tof - t)/dtdz

        n += 1

        if n >= 50:
            return None, None

    c = pred.stumpff_c(z)
    s = pred.stumpff_s(z)
    y = r1_norm + r2_norm - A * (1 - z * s) / np.sqrt(c)

    f = 1 - y/r1_norm
    g = A*np.sqrt(y/mu)
    gdot = 1 - y/r2_norm

    v1 = (r2 - f*r1)/g
    v2 = (gdot*r2 - r1)/g

    if np.cross(r1, v1)[2] < 0:
        if dir == 1:
            dir_n = 2
        else:
            dir_n = 1

        v1, v2 = lambert_uv(r1, r2, tof, mu, dir_n)

    return v1, v2


def lambert_mine(r1, r2, mu, tm=1):
    """
        :param r1: initial position
        :param r2: final position
        :param mu: gravitational parameter
        :param tm: transfer method (1 = short way, -1 = long way)

        :return amin: minimum semimajor axis
        :return emin: minimum eccentricity
        :return tmin_amin: minimum transfer time
        :return v1: initial velocity
    """
    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)
    cos_d_theta = np.dot(r1, r2)/(r1_norm*r2_norm)
    sin_d_theta = tm*np.sqrt(1 - cos_d_theta**2)
    c = np.sqrt(r1_norm**2 + r2_norm**2 - 2*r2_norm*r1_norm*cos_d_theta)
    s = (r1_norm + r2_norm + c)/2

    amin = s/2
    pmin = ((r1_norm*r2_norm)/c)*(1 - cos_d_theta)
    emin = np.sqrt(1 - (2*pmin)/s)

    alpha_e = np.pi
    beta_e = 2*np.arcsin(np.sqrt((s - c)/s))

    tmin_amin = np.sqrt((amin**3)/mu)*(alpha_e + tm*(beta_e - np.sin(beta_e)))
    v1 = np.sqrt(mu*pmin)/(r1_norm*r2_norm*sin_d_theta)*(r2 - (1 - (r2_norm/pmin)*(1 - cos_d_theta))*r1)

    return amin, emin, tmin_amin, v1


r1 = np.array([cns.re*1.2, 0, 0])
r2 = np.array([0, cns.re*2, 0])
amin, emin, tmin, v1 = lambert_mine(r1, r2, cns.mu, 1)

tspan = np.linspace(0.1, tmin, 1000)
rs = propagators.f_g_propagation(r1, v1, tspan, cns.mu)
plt.plot(rs[:, 0], rs[:, 1], 'k')

t2 = np.linspace(tmin, func.period(amin, cns.mu), 1000)
rs = propagators.f_g_propagation(r1, v1, t2, cns.mu)
plt.plot(rs[:, 0], rs[:, 1], 'k--')

x, y = aplt.plot_circle(cns.re, 100)
plt.fill(x, y, 'b', alpha=0.25)
plt.axis('equal')
plt.show()
