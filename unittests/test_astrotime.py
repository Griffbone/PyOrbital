import unittest
import astrotime as at


class TestFunctions(unittest.TestCase):
    def test_date_to_jd(self):
        tol = 0.00001   # equivalent to 0.864 seconds

        jdn, jdf = at.date_to_jd(2002, 3, 19, 12, 47, 0)
        ans = 2452353.03264
        self.assertTrue(ans - tol <= jdn + jdf <= ans + tol)

        jdn, jdf = at.date_to_jd(1996, 10, 26, 20, 14, 0)
        ans = 2450383.34306
        self.assertTrue(ans - tol <= jdn + jdf <= ans + tol)


if __name__ == '__main__':
    unittest.main()
