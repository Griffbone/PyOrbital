import numpy as np
from scipy.optimize import fsolve


def kepler_E(e, M):
    """Solve Kepler's equation for Eccentric anomaly
        :param e: eccentricity
        :param M: mean anomaly
        :return E: eccentric anomaly
    """

    if (e < 0) or (e > 1):
        print('ERROR: eccentricty must be below one')
        return None
    elif (M < 0) or (M > np.pi):
        print('ERROR: mean anomaly must be below one')
        return None

    return fsolve(lambda E: E - e*np.sin(E) - M, np.array([0]))

