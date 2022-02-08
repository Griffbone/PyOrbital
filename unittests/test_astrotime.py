import unittest
import astrotime as at


class TestFunctions(unittest.TestCase):
    def test_date_to_jd(self):
        tol = 0.00001

        JDN = at.date_to_jd(4, 2, 2022, 10, 37, 30)
        self.assertTrue(2459614.94271 - tol <= JDN <= 2459614.94271 + tol)

        JDN = at.date_to_jd(26, 10, 1996, 14, 20, 0)
        print(JDN)
        # self.assertTrue(2415021.02118 - tol <= JDN <= 2415021.02118 + tol)


if __name__ == '__main__':
    unittest.main()
