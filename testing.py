import propagators as prop
import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import cartopy.crs as crs
import astrotime as at
import functions as func
from tle import read_tle
from cartopy.feature.nightshade import Nightshade
from datetime import datetime, timezone, timedelta


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


rp = 200e3 + cns.re
ra = 485000e3 + cns.re
a = (rp + ra)/2
e = (ra - rp)/(ra + rp)


