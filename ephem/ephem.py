import astrotime as at
import functions as func
import constants as cns


class planetary_ephem:
    def __init__(self, a, e, i):
        pass

    def get_element(self, a0, a1, a2, a3, jd):
        """ Function to get an orbital element via polynomial interpolation
            :param a0: constant term
            :param a1: first coefficient
            :param a2: second coefficient
            :param a3: third coefficient
            :param jd: Julian date of epoch

            :return element: orbital element described by polynomial at epoch
        """

        T = at.j2000_cents(jd)
        element = a0 + a1*T + a2*T^2 + a3*T^3
        return element

mercury = planetary_ephem([0.387098310, 0, 0, 0], [0.20563175, 0.000020407, -0.0000000283, -0.00000000018])
