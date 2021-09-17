import numpy as np
import math
import time


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


print(plane_change(1000, 1000, 1000))
print(plane_change(1000, 1000))
