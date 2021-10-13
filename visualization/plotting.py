import matplotlib.pyplot as plt
import numpy as np
import pyOrbital.constants as cons


def perifocal_coords(a, e, thetas=np.linspace(0, 2*np.pi, 100), polar=True):
    """ Function to get polar perifocal coordinates of an orbit
        :param a: semimajor axis
        :param e: eccentricity
        :param thetas: array of true anaomalies
        :param polar: returns polar coordinates if true, cartesian coordinates if false
    """
    b = a*np.sqrt(1 - e**2)
    p = b**2/a

    rs = p/(1 - e*np.cos(thetas))

    if polar is True:
        return rs, thetas
    elif polar is False:
        return rs*np.cos(thetas), rs*np.sin(thetas)
    else:
        return None
        # raise TypeError("Expected <class 'bool'> for 'polar' argument but got {}".format(type(polar)))


def perifocal_plot(a, e, ax):
    p = (a*np.sqrt(1 - e**2))**2 / a

    thetas = np.linspace(0, np.pi*2, 1000)
    rs = p/(1 + e*np.cos(thetas))

    xs = rs*np.cos(thetas)
    ys = rs*np.sin(thetas)

    ax.plot(xs, ys)

# def perifocal_plot(a, e, plottype='cart', rc=cons.re, tas=np.array([])):
#     b = a*np.sqrt(1 - e**2)
#     p = b**2/a
#
#     thetas = np.linspace(0, np.pi*2, 100)
#     rs = p/(1 + e*np.cos(thetas))
#     rp = p/(1 + e*np.cos(tas))
#
#     if plottype.upper().replace(' ', '') in ['CART', 'CARTESIAN', 'XY']:
#         x = rs * np.cos(thetas)
#         y = rs * np.sin(thetas)
#         xp = rp * np.cos(tas)
#         yp = rp * np.sin(tas)
#
#         # fig, ax = plt.subplots(subplot_kw={'yeeter': 'True'})
#         plt.plot(rc * np.cos(thetas), rc * np.sin(thetas))
#         plt.plot(x, y)
#         plt.scatter(xp, yp)
#         plt.axis('equal')
#     elif plottype.upper().strip().replace(' ', '') in ['POLAR', 'RTHETA']:
#         fig, ax = plt.subplots(subplot_kw={'polar': 'True'})
#
#         ax.plot(thetas, np.ones(len(thetas))*rc)
#         ax.plot(thetas, rs)
#         ax.scatter(tas, rp)
#         ax.set_xticks(np.arange(0, np.pi*2, np.pi))
#         ax.set_yticks(np.linspace(0, ax.get_ylim()[1], 4, endpoint=False))
#         ax.set_rlabel_position(180)
#     else:
#         return None
#         # raise ValueError("unrecognized plot type: '{}'".format(plottype))

    # plt.title('Orbit Perifocal Plot; a={}, e={}'.format(a, e))
    # plt.show()


# perifocal_plot(cons.re*2, 0.3, plottype='xy', tas=np.linspace(0, np.pi/2, 10))
