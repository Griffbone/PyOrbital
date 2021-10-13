import numpy as np


def vis_viva(r1, r2, mu):
    """ Find velocity at apses of an orbit
        :param r1: first apsis radius
        :param r2: second apsis radius
        :param mu: parent body gravitational parameter

        :return v1: first apsis velocity
        :return v2: second apsis velocity
    """

    a = (r1 + r2)/2
    v1 = np.sqrt(mu*(2/r1 - 1/a))
    v2 = np.sqrt(mu*(2/r2 - 1/a))

    return v1, v2


def vis_viva_v(r, v, mu):
    """ Function to find the radius of the opposite apsis of an orbit from apsis radius and velocity
        :param r: radius of one apsis
        :param v: velocity of one apsis
        :param mu: parent body gravitational parameter

        :return r2: radius of opposite apsis
    """

    r2 = (-(v**2/mu - 2/r))**(-1)
    return r2


def plane_change(v1, di, v2=None):
    """ Function to calculate Delta-V for a plane change maneuver
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


def hohmann(r1, r2, mu):
    """ Function to perform a Hohmann transfer between two circular orbits
        :param r1: radius of first orbit
        :param r2: radius of second orbit
        :param mu: planet gravitational parameter
        :return dv1: delta-v of first burn
        :return dv2: delta-v of second burn
        :return T: period of the transfer orbit

        negative delta-v indicates the spacecraft must slow down (retrograde propulsion) to complete the transfer
    """

    v1 = np.sqrt(mu/r1)
    v2 = np.sqrt(mu/r2)

    vt1, vt2 = vis_viva(r1, r2, mu)

    dv1 = (vt1 - v1)
    dv2 = (v2 - vt2)

    T = 2*np.pi*np.sqrt(((r1+r2)/2)**3 / mu)

    return dv1, dv2, T


def impulsive_maneuver(ra, rp, dv, mu, pe=True):
    """ Function to calculate change in apses after an impulsive tangential maneuver
        :param ra: initial apoapsis radius
        :param rp: initial periapsis radius
        :param dv: burn delta-v
        :param mu: parent body gravitational parameter
        :param pe: true if boosting at periapsis
    """

    vp, va = vis_viva(ra, rp, mu)

    if pe is True :
        vp += dv
        ra = vis_viva_v(rp, vp, mu)
    else:
        va += dv
        rp = vis_viva_v(ra, va, mu)

    return ra, rp
