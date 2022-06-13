import astrotime as at
import functions as func
import constants as cns


class planetary_ephem:
    def __init__(self, elements, rates):
        self.elements = elements
        self.rates = rates

        pass

    def get_element(self, jd_et):
        """ Function to determine planetary orbital elements with linear rates
            :param jd_et: Julian date in ephemeris time

            :return elements: 6-element list of orbital elements (a, e, I, L, w bar, Omega)
        """
        T = at.j2000_cents(jd_et)
        elements = []

        for i in range(0, len(self.elements)):
            elements.append(self.elements[i] + self.rates[i]*T)




        return elements

        return element


# http://vadimchazov.narod.ru/text_pdf/XSChap8.pdf Table 8.10.2
mercury = planetary_ephem([0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593],
                          [0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081])
venus = planetary_ephem([0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255],
                        [0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418])
earth = planetary_ephem([1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0],
                        [0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0])
mars = planetary_ephem([1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891],
                       [0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343])