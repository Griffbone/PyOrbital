import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import constants as cons
import matplotlib.pyplot as plt

# STK propagators: https://help.agi.com/stk/index.htm#stk/vehSat_orbitProp_choose.htm
# IITB propagators: https://www.aero.iitb.ac.in/satelliteWiki/index.php/Orbit_Propagator


def kepler_ea(e, ma):
    """Solve Kepler's equation for Eccentric anomaly
        :param e: eccentricity
        :param ma: mean anomaly
        :return E: eccentric anomaly
    """

    if (e < 0) or (e > 1):
        print('ERROR: eccentricty must be below one')
        return None
    elif (ma < 0) or (ma > np.pi):
        print('ERROR: mean anomaly must be below one')
        return None

    return fsolve(lambda E: E - e*np.sin(E) - ma, np.array([0]))


def twobody_propagation(r, v, dt, n=1000):
    """ Propagate a satellite using two body propagation

    """

    def state_dot(time, state):
        x, y, z, x_dot, y_dot, z_dot = state

        r = np.array([x, y, z])
        rm = np.linalg.norm(r)

        gx = -((cons.mu*x)/rm)*(1 + 1.5*cons.j2)

        vx_dot = g[0]
        vy_dot = g[1]
        vz_dot = g[2]

        return np.array([x_dot, y_dot, z_dot, vx_dot, vy_dot, vz_dot])

    state_0 = np.array([r[0], r[1], r[2], v[0], v[1], v[2]])
    sol = solve_ivp(state_dot, [0, dt], state_0, 'DOP853', t_eval=np.linspace(0, dt, n))

    t = sol.t
    y = sol.y

    xs = y[0, :]
    ys = y[1, :]
    zs = y[2, :]

    return t, xs, ys, zs


t, x, y, z = twobody_propagation([cons.re, 0, 0], [0, 8000, 250], 90*60, n=10000)
plt.plot(x, y)
plt.axis('equal')
plt.show()


def j2_propagation(r, v, dt, n=1000):
    """ Propagate a satellite using two body propagation

    """

    def state_dot(time, state):
        x, y, z, x_dot, y_dot, z_dot = state

        r = np.array([x, y, z])
        g = -(cons.mu*r)/np.linalg.norm(r)**3

        vx_dot = g[0]
        vy_dot = g[1]
        vz_dot = g[2]

        return np.array([x_dot, y_dot, z_dot, vx_dot, vy_dot, vz_dot])

    state_0 = np.array([r[0], r[1], r[2], v[0], v[1], v[2]])
    sol = solve_ivp(state_dot, [0, dt], state_0, 'DOP853', t_eval=np.linspace(0, dt, n))

    t = sol.t
    y = sol.y

    xs = y[0, :]
    ys = y[1, :]
    zs = y[2, :]

    return t, xs, ys, zs