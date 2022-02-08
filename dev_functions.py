import numpy as np
import functions as func


def lla_to_ecef(lat, lon, alt):
    """ Function to convert geodetic coordinates and altitude to ECEF coordinates
        :param lat: geodetic latitude
        :param lon: geodetic longitude
        :param alt: altitude above ellipsoid (m)
    """
    e = 1/cns.flat
    lat = np.radians(lat)
    lon = np.radians(lon)

    xp = (cns.smae / np.sqrt(1 - e**2 * np.sin(lat)**2) + alt)*np.cos(lat)
    z = (cns.smae * (1 - e**2) / np.sqrt(1 - e**2 * np.sin(lat)**2) + alt)*np.sin(lat)

    x = xp*np.cos(lon)
    y = xp*np.sin(lon)

    return x, y, z


def lla_to_eci(lat, lon, alt, jdn):
    """ Function to convert geodetic coordinates and altitude to ECEF coordinates
        :param lat: geodetic latitude
        :param lon: geodetic longitude
        :param alt: altitude above ellipsoid (m)
        :param jdn: Julian day number
    """
    e = 1/cns.flat
    lat = np.radians(lat)
    lon = astt.theta_g(jdn) + np.radians(lon)

    xp = (cns.smae / np.sqrt(1 - e**2 * np.sin(lat)**2) + alt)*np.cos(lat)
    z = (cns.smae * (1 - e**2) / np.sqrt(1 - e**2 * np.sin(lat)**2) + alt)*np.sin(lat)

    x = xp*np.cos(lon)
    y = xp*np.sin(lon)

    return x, y, z


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
        Dt = func.t_between(sat.a, sat.e, sat.ta, ta_sat, sat.mu)

        ta_targ = func.ta_change(targ.a, targ.e, targ.ta, Dt, targ.mu)

        xs, ys, zs = sat.get_eci(ta_sat)
        xt, yt, zt = targ.get_eci(ta_targ)

        r = np.sqrt((xt - xs)**2 + (yt - ys)**2 + (zt - zs)**2)

        if r <= rs:
            intercept = True
            break

        if ta_sat > ta2:
            break

        ta_sat += dta

    return ta_sat, ta_targ, intercept


def patched_conics(sat, targ, rs, mu_t):
    ta_sat, ta_targ, incpt = patch_point(sat, targ, rs)

    if incpt is False:
        return None

    xs, ys, zs = func.elements_to_eci_pos(sat.a, sat.e, sat.i, sat.lan, sat.w, ta_sat)
    vxs, vys, vzs = func.elements_to_eci_vel(sat.a, sat.e, sat.i, sat.lan, sat.w, ta_sat)

    xt, yt, zt = func.elements_to_eci_pos(targ.a, targ.e, targ.i, targ.lan, targ.w, ta_targ)
    vxt, vyt, vzt = func.elements_to_eci_vel(targ.a, targ.e, targ.i, targ.lan, targ.w, ta_targ)

    vrel = np.array([vxs, vys, vzs]) - np.array([vxt, vyt, vzt])
    rrel = np.array([xs, ys, zs]) - np.array([xt, yt, zt])

    a, e, i, omega, w, M = func.vector_to_elements(rrel, vrel, mu_t)

    # print(M)
    return func.Elements(a, e, i, omega, w, 0, mu_t)
