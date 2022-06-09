import numpy as np
import constants as cns
import matplotlib.pyplot as plt
from plotting import plot_circle


# ====================== Helper Functions ======================
def stumpff_c(z):
    """ Evaluate the stumpff C series for given z
        :param z: z value
        :return c: c value
    """
    if z < 0:
        c = (1 - np.cosh(np.sqrt(-z)))/z
    else:
        c = (1 - np.cos(np.sqrt(z)))/z

    return c


def stumpff_s(z):
    """ Evaluate the Stumpff S series for given z
        :param z: z value
        :return s: s value
    """
    if z < 0:
        s = (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/np.sqrt((-z)**3)
    else:
        s = (np.sqrt(z) - np.sin(np.sqrt(z)))/np.sqrt(z**3)

    return s


def solve_uv(r, v, mu, tof, tol=1e-9):
    """ Function to solve for the universal variable chi from time of flight
        :param r: position vector
        :param v: velocity vector
        :param mu: gravitational parameter
        :param tof: time of flight
        :param tol: convergence tolerance

        :return xn: universal variable chi
    """
    r_norm = np.linalg.norm(r)
    sme = (1/2)*np.linalg.norm(v)**2 - mu/r_norm
    a = -mu/(2*sme)
    rvdot = np.dot(r, v)
    smu = np.sqrt(mu)

    xn = smu*tof/a
    zn = xn**2/a
    cn = stumpff_c(zn)
    sn = stumpff_s(zn)
    tn = (1/smu)*((rvdot/smu)*(xn**2)*cn + (1 - r_norm/a)*(xn**3)*sn + r_norm*xn)

    n = 0
    while abs(tn - tof) > tol:
        dtdx = (1/smu)*((xn**2*cn) + (rvdot/smu)*xn*(1 - zn*sn) + r_norm*(1 - zn*cn))

        xn = xn + (tof - tn)/dtdx
        zn = xn**2/a
        cn = stumpff_c(zn)
        sn = stumpff_s(zn)
        tn = (1/smu)*((rvdot/smu)*(xn**2)*cn + (1 - r_norm/a)*(xn**3)*sn + r_norm*xn)

        n += 1

    return xn


# ====================== Main Functions ======================
def f_g_state(r0, v0, mu, tof):
    """ Function to return final position and velocity using f and g functions
        :param r0: initial position
        :patam v0: initial velocity
        :param mu: gravitational parameter
        :param tof: time of flight

        :return r: final position
        :return v: final velocity
    """
    r0_norm = np.linalg.norm(r0)
    sme = (1/2)*np.linalg.norm(v0)**2 - mu/r0_norm
    a = -mu/(2*sme)

    x = solve_uv(r0, v0, mu, tof)
    z = x**2/a
    s = stumpff_s(z)
    c = stumpff_c(z)

    f = 1 - (x**2/r0_norm)*c
    g = tof - (x**3/np.sqrt(mu))*s
    r = f*r0 + g*v0
    r_norm = np.linalg.norm(r)
    gdot = 1 - (x**2/r_norm)*c
    fdot = np.sqrt(mu)/(r0_norm*r_norm)*x*(z*s - 1)
    v = fdot*r0 + gdot*v0

    return r, v
