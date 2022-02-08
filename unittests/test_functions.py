import unittest
import functions as func
import numpy as np
import constants as cns


class TestFunctions(unittest.TestCase):
    def test_vector_to_elements(self):
        c45 = np.cos(np.pi/4)
        mu = cns.mu/1000**3
        tol = 1e-6
        re = 6378
        vc = np.sqrt(mu/re)

        # Circular equatorial
        r = re*np.array([c45, c45, 0])
        v = vc*np.array([-c45, c45, 0])
        a, e, i, lan, w, ta, arglat, truelon, lonper = func.vector_to_elements(r, v, mu)

        self.assertTrue(6378 - tol <= a <= 6378 + tol)
        self.assertTrue(0 <= e <= tol)
        self.assertTrue(0 <= i <= tol)
        self.assertTrue(np.isnan(lan))
        self.assertTrue(np.isnan(w))
        self.assertTrue(np.isnan(ta))
        self.assertTrue(np.isnan(arglat))
        self.assertTrue(45 - tol <= truelon <= 45 + tol)
        self.assertTrue(np.isnan(lonper))

        # Retrograde circular equatorial
        r = re*np.array([-c45, c45, 0])
        v = vc*np.array([c45, c45, 0])
        a, e, i, lan, w, ta, arglat, truelon, lonper = func.vector_to_elements(r, v, mu)

        self.assertTrue(6378 - tol <= a <= 6378 + tol)
        self.assertTrue(0 <= e <= tol)
        self.assertTrue(180 - tol <= i <= 180 + tol)
        self.assertTrue(np.isnan(lan))
        self.assertTrue(np.isnan(w))
        self.assertTrue(np.isnan(ta))
        self.assertTrue(np.isnan(arglat))
        self.assertTrue(225 - tol <= truelon <= 225 + tol)
        self.assertTrue(np.isnan(lonper))

        # Retrograde circular at descending node
        r = re * np.array([0, -1, 0])
        v = vc * np.array([-c45, 0, -c45])
        a, e, i, lan, w, ta, arglat, truelon, lonper = func.vector_to_elements(r, v, mu)

        self.assertTrue(6378 - tol <= a <= 6378 + tol)
        self.assertTrue(0 <= e <= tol)
        self.assertTrue(135 - tol <= i <= 135 + tol)
        self.assertTrue(90 - tol <= lan <= 90 + tol)
        self.assertTrue(np.isnan(w))
        self.assertTrue(np.isnan(ta))
        self.assertTrue(180 - tol <= arglat <= 180 + tol)
        self.assertTrue(np.isnan(truelon))
        self.assertTrue(np.isnan(lonper))

        # Circular polar
        r = re * np.array([0, 0, 1])
        v = vc * np.array([c45, c45, 0])
        a, e, i, lan, w, ta, arglat, truelon, lonper = func.vector_to_elements(r, v, mu)

        self.assertTrue(6378 - tol <= a <= 6378 + tol)
        self.assertTrue(0 <= e <= tol)
        self.assertTrue(90 - tol <= i <= 90 + tol)
        self.assertTrue(225 - tol <= lan <= 225 + tol)
        self.assertTrue(np.isnan(w))
        self.assertTrue(np.isnan(ta))
        self.assertTrue(90 - tol <= arglat <= 90 + tol)
        self.assertTrue(np.isnan(truelon))
        self.assertTrue(np.isnan(lonper))

        # General circular
        r = re * np.array([-c45**2, c45**2, c45])
        v = vc * np.array([-c45, -c45, 0])
        a, e, i, lan, w, ta, arglat, truelon, lonper = func.vector_to_elements(r, v, mu)

        self.assertTrue(6378 - tol <= a <= 6378 + tol)
        self.assertTrue(0 <= e <= tol)
        self.assertTrue(45 - tol <= i <= 45 + tol)
        self.assertTrue(45 - tol <= lan <= 45 + tol)
        self.assertTrue(np.isnan(w))
        self.assertTrue(np.isnan(ta))
        self.assertTrue(90 - tol <= arglat <= 90 + tol)
        self.assertTrue(np.isnan(truelon))
        self.assertTrue(np.isnan(lonper))

        # Elliptical equatorial
        r = re * np.array([-c45, -c45, 0])
        v = 1.2 * vc * np.array([c45, -c45, 0])
        a, e, i, lan, w, ta, arglat, truelon, lonper = func.vector_to_elements(r, v, mu)

        self.assertTrue(11389.2857142 - tol <= a <= 11389.2857142 + tol)
        self.assertTrue(0.44 - tol <= e <= 0.44 + tol)
        self.assertTrue(0 - tol <= i <= 0 + tol)
        self.assertTrue(np.isnan(lan))
        self.assertTrue(np.isnan(w))
        self.assertTrue(0 - tol <= ta <= 0 + tol)
        self.assertTrue(np.isnan(arglat))
        self.assertTrue(np.isnan(truelon))
        self.assertTrue(225 - tol <= lonper <= 225 + tol)

        # General elliptical
        r = re * np.array([0, c45, c45])
        v = 1.2 * vc * np.array([-1, 0, 0])
        a, e, i, lan, w, ta, arglat, truelon, lonper = func.vector_to_elements(r, v, mu)

        self.assertTrue(11389.2857142 - tol <= a <= 11389.2857142 + tol)
        self.assertTrue(0.44 - tol <= e <= 0.44 + tol)
        self.assertTrue(45 - tol <= i <= 45 + tol)
        self.assertTrue(0 - tol <= lan <= 0 + tol)
        self.assertTrue(90 - tol <= w <= 90 + tol)
        self.assertTrue(0 - tol <= ta <= 0 + tol)
        self.assertTrue(90 - tol <= arglat <= 90 + tol)
        self.assertTrue(np.isnan(truelon))
        self.assertTrue(np.isnan(lonper))

        # General parabolic
        r = re * np.array([0, c45, c45])
        v = np.sqrt(2) * vc * np.array([-1, 0, 0])
        a, e, i, lan, w, ta, arglat, truelon, lonper = func.vector_to_elements(r, v, mu)

        self.assertTrue(a == np.inf)
        self.assertTrue(1 - tol <= e <= 1 + tol)
        self.assertTrue(45 - tol <= i <= 45 + tol)
        self.assertTrue(0 - tol <= lan <= 0 + tol)
        self.assertTrue(90 - tol <= w <= 90 + tol)
        self.assertTrue(0 - tol <= ta <= 0 + tol)
        self.assertTrue(90 - tol <= arglat <= 90 + tol)
        self.assertTrue(np.isnan(truelon))
        self.assertTrue(np.isnan(lonper))

        # Radial orbit
        r = re * np.array([0, c45, c45])
        v = 1.2 * vc * np.array([0, c45, c45])
        a, e, i, lan, w, ta, arglat, truelon, lonper = func.vector_to_elements(r, v, mu)

        self.assertTrue(np.isnan(a))
        self.assertTrue(np.isnan(e))
        self.assertTrue(np.isnan(i))
        self.assertTrue(np.isnan(lan))
        self.assertTrue(np.isnan(w))
        self.assertTrue(np.isnan(ta))
        self.assertTrue(np.isnan(arglat))
        self.assertTrue(np.isnan(truelon))
        self.assertTrue(np.isnan(lonper))

        pass

    def test_elements_to_vector(self):
        c45 = np.cos(np.pi/4)
        mu = cns.mu
        tol = 1e-6

        # Circular equatorial
        rc, vc = func.elements_to_vector(6378e3, 0, 0, 0, 0, 45, mu)
        r = 6378e3*np.array([c45, c45, 0])
        v = np.sqrt(mu/6378e3)*np.array([-c45, c45, 0])

        self.assertTrue((r - tol <= rc).all() and (rc <= r + tol).all())
        self.assertTrue((v - tol <= vc).all() and (vc <= v + tol).all())

        # Retrograde circular equatorial
        rc, vc = func.elements_to_vector(6378e3, 0, 180, 0, 0, 135, mu)
        r = 6378e3*np.array([-c45, -c45, 0])
        v = np.sqrt(mu/6378e3)*np.array([-c45, c45, 0])

        self.assertTrue((r - tol <= rc).all() and (rc <= r + tol).all())
        self.assertTrue((v - tol <= vc).all() and (vc <= v + tol).all())

        # Retrograde circular at descending node
        # rc, vc = func.elements_to_vector(6378e3, 0, 135, 90, 0, 180, mu)
        # r = 6378e3*np.array([0, -1, 0])
        # v = np.sqrt(mu/6378e3)*np.array([-c45, 0, -c45])
        #
        # self.assertTrue((r - tol <= rc).all() and (rc <= r + tol).all())
        # self.assertTrue((v - tol <= vc).all() and (vc <= v + tol).all())

        # Circular polar
        # rc, vc = func.elements_to_vector(6378e3, 0, 90, 45, 0, 0, mu)
        # r = 6378e3*np.array([c45, c45, 0])
        # v = np.sqrt(mu/6378e3)*np.array([0, 0, 1])
        #
        # self.assertTrue((r - tol <= rc).all() and (rc <= r + tol).all())
        # self.assertTrue((v - tol <= vc).all() and (vc <= v + tol).all())

        # General circular
        # rc, vc = func.elements_to_vector(6378e3, 0, 45, 45, 0, 90, mu)
        # r = 6378e3*np.array([-c45**2, c45**2, c45])
        # v = np.sqrt(mu/6378e3)*np.array([-c45, -c45, 0])
        #
        # self.assertTrue((r - tol <= rc).all() and (rc <= r + tol).all())
        # self.assertTrue((v - tol <= vc).all() and (vc <= v + tol).all())

        #
        # # general elliptical (Vallado ex. 2-5, using lower tolerances here)
        # r = 1000*np.array([6524.834, 6862.875, 6448.296])
        # v = 1000*np.array([4.901327, 5.533756, -1.976341])
        # a, e, i, lan, w, ta = func.vector_to_elements(r, v, mu)
        # tol = 10
        # self.assertTrue(36127.343e3 - tol <= a <= 36127.343e3 + tol)
        # tol = 0.1
        # self.assertTrue(0.832853 - tol <= e <= 0.832853 + tol)
        # self.assertTrue(87.870 - tol <= i <= 87.870 + tol)
        # self.assertTrue(227.898 - tol <= lan <= 227.898 + tol)
        # self.assertTrue(53.38 - tol <= w <= 53.38 + tol)
        # self.assertTrue(92.335 - tol <= ta <= 92.335 + tol)
        # tol = 1e-6
        #
        # # elliptical equatorial
        #
        # # parabolic
        # r = 6378e3*np.array([c45, c45, 0])
        # v = np.sqrt(2*mu/6378e3)*np.array([-c45, c45, 0])
        # a, e, i, lan, w, ta = func.vector_to_elements(r, v, mu)
        # self.assertTrue(a == np.inf)
        # self.assertTrue(1 - tol <= e <= 1 + tol)
        # self.assertTrue(0 <= i <= tol)
        # self.assertTrue(0 <= lan <= tol)
        # self.assertTrue(45 - tol <= w <= 45 + tol)
        # self.assertTrue(0 <= ta <= tol)
        #
        # # hyperbolic



if __name__ == '__main__':
    unittest.main()
