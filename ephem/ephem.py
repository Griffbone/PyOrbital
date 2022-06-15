import numpy

import astrotime as at
import functions as func
import constants as cns
import numpy as np
import prediction as pred
import matplotlib.pyplot as plt
import gauss
import propagators
import plotting as aplt


class planetary_ephem:
    def __init__(self, elements, rates, color='k'):
        self.elements = elements
        self.rates = rates
        self.color = color

        pass

    def get_elements(self, jd_et):
        """ Function to determine planetary orbital elements with linear rates
            :param jd_et: Julian date in ephemeris time

            :return elements: 6-element list of orbital elements (a, e, I, L, w bar, Omega)
        """
        T = at.j2000_cents(jd_et)

        a = (self.elements[0] + self.rates[0]*T)*cns.au
        e = (self.elements[1] + self.rates[1]*T)
        I = func.wrap_to_360(self.elements[2] + self.rates[2]*T)
        L = func.wrap_to_360(self.elements[3] + self.rates[3]*T)
        wbar = func.wrap_to_360(self.elements[4] + self.rates[4]*T)
        Omega = func.wrap_to_360(self.elements[5] + self.rates[5]*T)

        return a, e, I, L, wbar, Omega

    def get_state(self, jd_et):
        """ Function to determine body state in J2000 frame
            :param jd_et: Ephemeris time Julian date

            :return r_j2000: body position in J2000 mean ecliptic/equinox frame (m)
            :return v_j2000: body velocity in J2000 mean ecliptic/equinox frame (m/s)
        """
        a, e, I, L, wbar, Omega = self.get_elements(jd_et)
        w = wbar - Omega
        M = np.deg2rad((L - wbar) % 360)
        E = func.kepler_inverse(M, e)

        xp = a*(np.cos(E) - e)
        yp = a*np.sqrt(1 - e**2)*np.sin(E)
        r_j2000 = func.perifocal_to_eci(xp, yp, np.deg2rad(Omega), np.deg2rad(I), np.deg2rad(w))

        return r_j2000

    def get_state_array(self, jd_ets):
        """ Function to get array of body states over given time span
            :param jd_ets: 1D array of ephemeris time julian dates

            :return rs: body positions in heliocentric J2000 frame (m)
        """
        rs = np.empty((0, 3))

        for jd in jd_ets:
            r = self.get_state(jd)
            rs = np.vstack([rs, r])

        return rs


# http://vadimchazov.narod.ru/text_pdf/XSChap8.pdf Table 8.10.2
# https://ssd.jpl.nasa.gov/planets/approx_pos.html
mercury = planetary_ephem([0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593],
                          [0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081], 'grey')
venus = planetary_ephem([0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255],
                        [0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418], 'y')
earth = planetary_ephem([1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0],
                        [0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0], 'b')
mars = planetary_ephem([1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891],
                       [0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343], 'r')
jupiter = planetary_ephem([5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909],
                          [-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106], 'g')
saturn = planetary_ephem([9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448],
                         [-0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794], 'y')
uranus = planetary_ephem([19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503],
                         [-0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589], 'b')
neptune = planetary_ephem([30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574],
                          [0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664], 'b')


def viz_system(jd_ets, ephems, au=False):
    """ Function to visualize the solar system
        :param jd_ets: 1D array of ephemeris time Julian dates
        :param ephems: ephemeris objects to visualize
        :param au: plot values in AU or m
    """

    if au is True:
        coeff = 1/cns.au
    else:
        coeff = 1

    for body in ephems:
        rs = body.get_state_array(jd_ets)*coeff
        plt.plot(rs[:, 0], rs[:, 1], color=body.color)
        plt.scatter([rs[0][0]], [rs[0][1]], marker='x', color=body.color)
        plt.scatter([rs[-1][0]], [rs[-1][1]], marker='o', color=body.color)

    plt.scatter([0], [0], color='y', marker='*')
    plt.axis('equal')


def transfer(src, dest, jd0, jdf):
    """ Function to determine an orbit transfer between celestial objects
        :param src: source body
        :param dest: destination body
        :param jd0: julian date of departure
        :param jdf: julian date of arrival

        :return vinf_1: hyperbolic excess speed at departure
        :return vinf_2: hyperbolic excess speed at arrival
    """
    rs, vs = src.get_pos(jd0)
    rd, vd = dest.get_pos(jdf)

    v1, v2 = gauss.lambert_uv(rs, rd, (jdf - jd0)*24*60*60, cns.mu_sun)

    if v1 is None or v2 is None:
        vinf_1 = np.nan
        vinf_2 = np.nan
    else:
        vinf_1 = np.linalg.norm(v1 - vs)
        vinf_2 = np.linalg.norm(v2 - vs)

    # ts = np.linspace(0, (jdf - jd0)*24*60**2)
    # xs = []
    # ys = []
    # for t in ts:
    #     r, _ = pred.f_g_state(rs, v1, cns.mu_sun, t)
    #     xs.append(r[0])
    #     ys.append(r[1])

    # plt.plot(xs, ys, 'k')

    return vinf_1, vinf_2


def porkchop(src, dest, jdmin, jdmax, amin, amax, res=100):
    """ Function to create data for a porkchop plot of transfer opportunities between two bodies
        :param src: source body
        :param dest: destination body
        :param jdmin: min julian date
        :param jdmax: max julian date
        :param amin: minimum arrival julian date
        :param amax: maximum arrival julian date

        :return xx: x-axis meshgrid of julian dates
        :return yy: y-axis meshgrid of time of flights
        :return zz: z-axis meshgrid of departure C3s
    """
    x = np.linspace(jdmin, jdmax, res)
    y = np.linspace(amin, amax, res)
    xx, yy = np.meshgrid(x, y)
    zz = np.empty((len(y), len(x)))

    n = 0

    for i in range(0, len(y)):
        for j in range(0, len(x)):
            vinf_1, _ = transfer(src, dest, x[j], y[i])
            zz[i, j] = (vinf_1/1000)**2

            n += 1
            print(n/(len(x)*len(y)))

    return xx, yy, zz


def viz_transfer(r1, v1, tof):
    ts = np.linspace(0, tof, 1000)
    xs = []
    ys = []

    for t in ts:
        r, _ = pred.f_g_state(r1, v1, cns.mu_sun, t)
        xs.append(r[0])
        ys.append(r[1])

    plt.plot(xs, ys, 'k')


jd0, _ = at.date_to_jd(1998, 11, 29, 0, 0, 0)
# jdf = jd0 + 365
jdf, _ = at.date_to_jd(2022, 6, 15, 0, 0, 0)
tspan = np.linspace(jd0, jdf, 10000)
print(at.jd_to_day(jdf))

viz_system(tspan, [mercury, venus, earth, mars, jupiter, saturn, uranus, neptune], True)
plt.show()
