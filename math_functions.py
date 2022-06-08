"""
    Functions for basic mathematics.
    Author      : Griffin Jourda
    Date        : 2/16/22
    Last Edited : 2/16/22

    Functions
        wrap_to_360 :   wrap an angle in degrees to range [0, 360]
        wrap_to_2pi :   wrap an angle in radians to range [0, 2*pi]

    Dependencies
        numpy       :   pi
"""

import numpy as np


def wrap_to_360(theta):
    """ Wrap an angle in degrees into range range [0, 360] deg

        :param theta: angle to wrap (deg)
        :return theta_w: wrapped angle (deg)
    """

    theta_w = theta % 360
    return theta_w


def wrap_to_2pi(theta):
    """ Wrap an angle in radians into range [0, 2*pi]

        :param theta: angle to wrap (rad)
        :return theta_w: wrapped angle (rad)
    """

    theta_w = theta % (2*np.pi)
    return theta_w
