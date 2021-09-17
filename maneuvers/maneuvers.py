import numpy as np
import math
import time
import pyOrbital.constants as cons


def plane_change(v1, di, v2=None):
    """ Function to calculate Delta-V for a simple plane change maneuver
        :param v1: initial vehicle velocity
        :param di: inclination change or velocity vector direction change
        :param v2: final vehicle velocity
        :return dv: delta-v of the maneuver
    """

    if v2 is None:
        dv = 2*v1*np.sin(di/2)
    else:
        dv = np.sqrt(2*v1**2*(1 - np.cos(di)))

    return abs(dv)


def vis_viva(ra, rp, mu):
    """ Find velocity at apses of an orbit
        :param ra: apoapsis velocity
        :param rp: periapsis velocity
        :param mu: planet gravitational parameter

        :return va: apoapsis velocity
        :return vp: periapsis velocity
    """

    a = (ra + rp)/2
    va = np.sqrt(mu*(2/ra - 1/a))
    vp = np.sqrt(mu*(2/rp - 1/a))

    return va, vp


def hohmann(r1, r2, mu):
    """ Function to perform a Hohmann transfer between two circular orbits
        :param r1: radius of first orbit
        :param r2: radius of second orbit
        :param mu: planet gravitational parameter
        :return dv1: delta-v of first burn
        :return dv2: delta-v of second burn

        negative delta-v indicates the spacecraft must slow down (retrograde propulsion) to complete the transfer
    """

    vi = np.sqrt(mu/r1)
    vf = np.sqrt(mu/r2)
    a = (r1 + r2)/2
    vtp = np.sqrt(mu*(2/r1 - 1/a))
    vta = np.sqrt(mu*(2/r2 - 1/a))

    dv1 = vtp - vi
    dv2 = vf - vta

    return dv1, dv2

# print(vis_viva(cons.re + 100000, cons.re + 500000, cons.mu))

print(hohmann(590000 + cons.re, 650000 + cons.re, cons.mu))

