import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import constants as cons
import matplotlib.pyplot as plt
import functions as func

# STK propagators: https://help.agi.com/stk/index.htm#stk/vehSat_orbitProp_choose.htm
# IITB propagators: https://www.aero.iitb.ac.in/satelliteWiki/index.php/Orbit_Propagator


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


def kepler_propagation(a, e, i, lan, w, ta, dt, n=1000):
    """ Function to propagate an orbit using Kepler's equation
        :param a: semimajor axis (m)
        :param e: eccentricity
        :param i: inclination  (deg)
        :param lan: longitude of ascending node (deg)
        :param w: argument of periapsis (deg)
        :param ta: true anomaly
        :param dt: final time (s)

    """

    ts = np.linspace(0, dt, n)
    n = np.sqrt(cons.mu/a**3)
    mas = n*ts

    tas = []
    xs = []
    ys = []
    zs = []

    for ma in mas:
        ta = kepler_ta(e, ma)

    return ts, np.array(tas)


def twobody_propagation(r, v, dt, n=1000, step=1e3):
    """ Propagate a satellite using two body propagation
        :param r: initial position (m)
        :param v: initial velocity (m/s)
        :param dt: end time of simulation (s)
        :param n: number of data points to return
        :param step: maximum step size (s)
    """

    def state_dot(time, state):
        x, y, z, x_dot, y_dot, z_dot = state

        r = np.array([x, y, z])
        rm = np.linalg.norm(r)

        g = -(cons.mu*r)/rm**3

        vx_dot = g[0]
        vy_dot = g[1]
        vz_dot = g[2]

        return np.array([x_dot, y_dot, z_dot, vx_dot, vy_dot, vz_dot])

    state_0 = np.array([r[0], r[1], r[2], v[0], v[1], v[2]])
    sol = solve_ivp(state_dot, [0, dt], state_0, 'DOP853', t_eval=np.linspace(0, dt, n), max_step=step)

    t = sol.t
    y = sol.y

    xs = y[0, :]
    ys = y[1, :]
    zs = y[2, :]

    return t, xs, ys, zs


def j2_propagation(r, v, dt, n=1000, step=1e3):
    """ Propagate a satellite using two body propagation with J2 perturbation
        :param r: initial position (m)
        :param v: initial velocity (m/s)
        :param dt: end time of simulation (s)
        :param n: number of data points to return
        :param step: maximum step size (s)
    """

    def state_dot(time, state):
        x, y, z, x_dot, y_dot, z_dot = state

        r = np.sqrt(x**2 + y**2 + z**2)

        gx = -(cons.mu*x/(r**3))*(1 + 1.5*cons.j2*(cons.re**2/r**2) - 7.5*cons.j2*(cons.re**2*z**2)/r**4)
        gy = -(cons.mu*y/(r**3))*(1 + 1.5*cons.j2*(cons.re**2/r**2) - 7.5*cons.j2*(cons.re**2*z**2)/r**4)
        gz = -(cons.mu*z/(r**3))*(1 + 4.5*cons.j2*(cons.re**2/r**2) - 7.5*cons.j2*(cons.re**2*z**2)/r**4)

        return np.array([x_dot, y_dot, z_dot, gx, gy, gz])

    state_0 = np.array([r[0], r[1], r[2], v[0], v[1], v[2]])
    sol = solve_ivp(state_dot, [0, dt], state_0, 'DOP853', t_eval=np.linspace(0, dt, n), max_step=step)

    t = sol.t
    y = sol.y

    xs = y[0, :]
    ys = y[1, :]
    zs = y[2, :]
    vxs = y[3, :]
    vys = y[4, :]
    vzs = y[5, :]

    return t, xs, ys, zs, vxs, vys, vzs
