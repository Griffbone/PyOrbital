import numpy

import astrotime as at
import functions as func
import constants as cns
import numpy as np
import prediction as pred
import matplotlib.pyplot as plt
import gauss
import propagators


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
        i = func.wrap_to_360(self.elements[2] + self.rates[2]*T)
        L = func.wrap_to_360(self.elements[3] + self.rates[3]*T)
        wbar = func.wrap_to_360(self.elements[4] + self.rates[4]*T)
        Omega = func.wrap_to_360(self.elements[5] + self.rates[5]*T)

        return a, e, i, L, wbar, Omega

    def get_coes(self, jd_et):
        a, e, i, L, wbar, Omega = self.get_elements(jd_et)

        w = wbar - Omega
        M = L - wbar
        ta = func.kepler_ta(e, M)

        return a, e, i, Omega, w, ta

    def get_pos(self, jd_et):
        a, e, i, Omega, w, ta = self.get_coes(jd_et)
        r, v = func.elements_to_vector(a, e, i, Omega, w, ta, None, None, None, cns.mu_sun)

        return r, v

    def get_state_array(self, jd_ets):
        rs = np.empty((0, 3))
        vs = np.empty((0, 3))

        for jd in jd_ets:
            a, e, i, Omega, w, ta = self.get_coes(jd)
            r, v = func.elements_to_vector(a, e, i, Omega, w, ta, None, None, None, cns.mu_sun)
            rs = np.vstack([rs, r])
            vs = np.vstack([vs, v])

        return rs, vs


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


def viz_system(jds, ephems):
    """ Function to visualize the solar system
    """

    for body in ephems:
        rs, _ = body.get_state_array(jds)
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


jd0, _ = at.date_to_jd(2005, 6, 1, 0, 0, 0)
jdf, _ = at.date_to_jd(2005, 10, 1, 0, 0, 0)

amin, _ = at.date_to_jd(2006, 1, 1, 0, 0, 0)
amax, _ = at.date_to_jd(2006, 11, 1, 0, 0, 0)

x, y, z = porkchop(earth, mars, jd0, jdf, amin, amax, 100)
# print(z)

plt.contour(x, y, z, levels=np.linspace(0, 40, 20))
plt.show()
