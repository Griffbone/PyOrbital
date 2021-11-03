import numpy as np
import functions as func


def ground_track(t0, ts, xs, ys, zs):
    """ Function to return geodetic lat/lon from ECI positions
        :param t0: epoch (JDN)
        :param ts: array of times since epoch (s)
        :param xs: ECI X position
        :param ys: ECI Y position
        :param zs: ECI Z position
    """

    phis = []
    lams = []

    for i in range(0, len(ts)):
        t = ts[i]/60**2
        phi, lam = func.eci_to_ll(xs[i], ys[i], zs[i], t0, t)

        phis.append(phi)
        lams.append(lam)

    return phis, lams


def patch_point(sat, targ, rs):
    """ Function to determine the patch point for a patched-conic problem
        :param sat_elements: satelltie orbital elements
        :param targ_elements: target orbital elements
        :param rs: target SOI radius
    """

    if targ.mu != sat.mu:
        raise ValueError('Target and satellite do not orbit the same body')

    # Determine if ap is high enough for intersection
    if sat.ap < targ.pe - rs:
        return None

    # Determine intersection TA pair
    ta1 = np.degrees(np.arccos(sat.a*(1 - sat.e**2)/((targ.a - rs)*sat.e) - 1/sat.e))
    ta2 = 360 - ta1

    # Step between these TAs to find intersection
    dta = (ta2 - ta1)/1000
    ta_sat = ta1

    intercept = False

    while True:
        Dt = func.t_between(sat.a, sat.e, ta1, ta_sat, sat.mu)

        ta_targ = func.ta_change(targ.a, targ.e, targ.ta, Dt, targ.mu)

        xs, ys, zs, _, _, _ = sat.get_eci(ta_sat)
        xt, yt, zt, _, _, _ = targ.get_eci(ta_targ)

        r = np.sqrt((xt - xs)**2 + (yt - ys)**2 + (zt - zs)**2)

        if r <= rs:
            intercept = True
            break

        if ta_sat > ta2:
            break

        ta_sat += dta

    return ta_sat, ta_targ, intercept


def