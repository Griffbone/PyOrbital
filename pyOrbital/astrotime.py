from constants import j2000
import numpy as np


def date_to_jd(day, month, year, hour=0, mins=0, sec=0):
    """ Function to calculate Julian Day from Gregorian Calendar date
        Follows Curtis Eq. 5.48

        :param day: gregorian calendar day
        :param month: gregorian calendar month
        :param year: gregorian calendar year
        :param hour: hour of day
        :param mins: minute of hour
        :param sec: second of minute
        :return jd: Julian date (days)
    """

    if year < 1901 or year > 2099:
        print('Cannot calculate JD for years before 1901 or after 2099')
        return None
    elif month < 1 or month > 12:
        print('Cannot calculate JD for months less than 1 or greater than 12')
        return None
    elif day < 1 or day > 31:
        print('Cannot calculate JD for days less than 1 or greater than 31')
        return None

    jd = 367*year - int(7*(year + int((month + 9)/12))/4) + int((275*month)/9) + day + 1721013.5
    ut = hour + mins/60 + sec/3600
    jd += ut/24

    return jd


def jd_to_date(jd):
    """ Function to calculate Gregorian calendar date from Julian date

        :param jd: julian date (days)
        :return date: gregorian calendar date
        :type date: tuple (day, month, year, hour, min, sec)
    """

    pass


def jd_to_t0(jd, a=0):
    """ Function to return Julian centuries since J2000 epoch
        Follows Curtis Eqs. 5.49 - 5.52

        :param jd: Julian date (days)
        :param a: east longitude (deg)
        :return T0: Julian centuries since J2000
        :return theta_g: Greenwich sidereal time (deg)
        :return theta: local sidereal time (deg)
    """

    jd2 = np.floor(jd) + 0.5
    ut = jd - jd2
    jd = jd2

    T0 = (jd - j2000)/36525
    theta_g = 100.4606184 + 36000.77004*T0 + 0.000387933*T0**2 - 2.583e-8*T0**3

    theta_g = theta_g % 360 + 360.98564724*ut
    theta = theta_g + a
    theta = theta % 360

    return T0, theta_g, theta
