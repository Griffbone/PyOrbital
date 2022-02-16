import numpy as np
from scipy.optimize import fsolve

mercury = {
    'name': 'mercury',

    # semimajor axis
    'a': 0.38709927,
    'adot': 0.00000037,

    # eccentricity
    'e': 0.20563593,
    'edot': 0.00001906,

    # inclination
    'I': 7.00497902,
    'Idot': -0.00594749,

    # mean longitude
    'L': 252.25032350,
    'Ldot': 149472.67411175,

    # longitude of perihelion
    'w': 77.45779628,
    'wdot': 0.16047689,

    # longitude of ascending node
    'Lan': 48.33076593,
    'Landot': -0.12534081,

    # positions:
    'x': [],
    'y': [],
    'z': []
}

venus = {
    'name': 'venus',

    # semimajor axis
    'a': 0.72333566,
    'adot': 0.00000390,

    # eccentricity
    'e': 0.00677672,
    'edot': -0.00004107,

    # inclination
    'I': 3.39467605,
    'Idot': -0.00078890,

    # mean longitude
    'L': 181.97909950,
    'Ldot': 58517.81538729,

    # longitude of perihelion
    'w': 131.60246718,
    'wdot': 0.00268329,

    # longitude of ascending node
    'Lan': 76.67984255,
    'Landot': -0.27769418,

    # positions:
    'x': [],
    'y': [],
    'z': []
}

earth = {
    'name': 'earth',

    # semimajor axis
    'a': 1.00000261,
    'adot': 0.00000562,

    # eccentricity
    'e': 0.01671123,
    'edot': -0.00004392,

    # inclination
    'I': -0.00001531,
    'Idot': -0.01294668,

    # mean longitude
    'L': 100.46457166,
    'Ldot': 35999.37244981,

    # longitude of perihelion
    'w': 102.93768193,
    'wdot': 0.32327364,

    # longitude of ascending node
    'Lan': 0,
    'Landot': 0.0,

    # positions:
    'x': [],
    'y': [],
    'z': []
}

mars = {
    'name': 'mars',

    # semimajor axis
    'a': 1.52371034,
    'adot': 0.00001847,

    # eccentricity
    'e': 0.09339410,
    'edot': 0.00007882,

    # inclination
    'I': 1.84969142,
    'Idot': -0.00813131,

    # mean longitude
    'L': -4.55343205,
    'Ldot': 19140.30268499,

    # longitude of perihelion
    'w': -23.94362959,
    'wdot': 0.44441088,

    # longitude of ascending node
    'Lan': 49.55953891,
    'Landot': -0.29257343,

    # positions:
    'x': [],
    'y': [],
    'z': []
}

jupiter = {
    'name': 'Jupiter',

    # semimajor axis
    'a': 5.20288700,
    'adot': -0.00011607,

    # eccentricity
    'e': 0.04838624,
    'edot': -0.00013253,

    # inclination
    'I': 1.30439695,
    'Idot': -0.00183714,

    # mean longitude
    'L': 34.39644051,
    'Ldot': 3034.74612775,

    # longitude of perihelion
    'w': 14.72847983,
    'wdot': 0.21252668,

    # longitude of ascending node
    'Lan': 100.47390909,
    'Landot': 0.20469106,

    # positions:
    'x': [],
    'y': [],
    'z': []
}

pluto = {
    'name': 'Pluto',

    # semimajor axis
    'a': 39.48211675,
    'adot': -0.0031596,

    # eccentricity
    'e': 0.24882730,
    'edot': 0.00005170,

    # inclination
    'I': 17.14001206,
    'Idot': 0.00004818,

    # mean longitude
    'L': 238.92903833,
    'Ldot': 145.20780515,

    # longitude of perihelion
    'w': 224.06891629,
    'wdot': -0.04062942,

    # longitude of ascending node
    'Lan': 110.30393684,
    'Landot': -0.01183482,

    # positions:
    'x': [],
    'y': [],
    'z': []
}


def get_elements_from_ephemeris(body):
    a = body['a']
    e = body['e']
    i = body['I']
    w = body['w'] - body['Lan']
    lan = body['Lan']

    print([w, lan])
    M = body['L'] - body['w']

    return a, e, i, w, lan, M


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


def norm_angle_deg(theta):
    return (theta + 180) % 360 - 180


def ephemeris(planet, T):
    T = (T - 2451545.0)/36525   # time past J2000

    a = planet['a'] + planet['adot']*T
    e = planet['e'] + planet['edot']*T
    I = planet['I'] + planet['Idot']*T
    L = planet['L'] + planet['Ldot']*T
    wbar = planet['w'] + planet['wdot']*T
    Lan = planet['Lan'] + planet['Landot']*T

    w = wbar - Lan          # argument of perihelion
    M = L - wbar            # mean anomaly
    M = norm_angle_deg(M)   # mean anomaly wrapped to [-180, 180]
    estar = (180/np.pi)*e   # eccentricity

    # solve Kepler's equation for eccentric anomaly
    E0 = M + estar*np.sin(np.radians(M))
    E = E0
    dE = 1

    while dE >= 1e-6:
        dM = M - (E - estar*np.sin(np.radians(E)))
        dE = dM/(1 - e*np.cos(np.radians(E)))
        E += dE

    # heliocentric coordinates
    xp = a*(np.cos(np.radians(E)) - e)
    yp = a*np.sqrt(1 - e**2)*np.sin(np.radians(E))
    zp = 0

    w = np.radians(w)
    Lan = np.radians(Lan)
    I = np.radians(I)

    recl = np.array([
        (np.cos(w)*np.cos(Lan) - np.sin(w)*np.sin(Lan)*np.cos(I))*xp +
        (-np.sin(w)*np.cos(Lan) - np.cos(w)*np.sin(Lan)*np.cos(I))*yp,

        (np.cos(w)*np.sin(Lan) + np.sin(w)*np.cos(Lan)*np.cos(I))*xp +
        (-np.sin(w)*np.sin(Lan) + np.cos(w)*np.cos(Lan)*np.cos(I))*yp,

        np.sin(w)*np.sin(I)*xp + np.cos(w)*np.sin(I)*yp
    ])

    return recl


# planets = [mercury, venus, earth, mars]
# n = 1
# Ts = np.linspace(2451545, 2451545 + 365*n, int(n*100))
#
# for T in Ts:
#     for body in planets:
#         r = ephemeris(body, T)
#         body['x'].append(r[0])
#         body['y'].append(r[1])
#         body['z'].append(r[2])
#
# fig, ax = plt.subplots(1, 3)
# leg = []
# for p in planets:
#     ax[0].plot(p['x'], p['y'])
#     ax[1].plot(p['y'], p['z'])
#     ax[2].plot(p['x'], p['z'])
#
#     leg.append(p['name'])
#
# leg.append('Sun')
# ax[0].scatter(0, 0, 0, color='y')
# ax[0].set_aspect('equal')
#
# ax[1].scatter(0, 0, 0, color='y')
# ax[1].set_aspect('equal')
#
# ax[2].scatter(0, 0, 0, color='y')
# ax[2].set_aspect('equal')
#
# fig.legend(leg)
# plt.show()
