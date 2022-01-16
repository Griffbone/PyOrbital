import functions
from constants import j2000
import datetime


def date_to_jd(day, month, year, h, m, s):
    """ Function to calculate Julian Day from Gregorian Calendar date
        Follows Curtis Eq. 5.48

        :param day: gregorian calendar day
        :param month: gregorian calendar month
        :param year: gregorian calendar year
        :return jd: Julian date (days)
    """

    jd = 367*year - int((7/4) * (year + int((month + 9)/12))) + int((275/9)*month) + day + 1721013.5
    ut =  h/24 + m/(24*60) + s/(24*60**2)

    return jd + ut


def datetime_to_jd(date):
    """ Function to calculate Julian Day from datetime object
        Follows Curtis Eq. 5.48

        :param date: datetime object
        :return jd: Julian date (days)
    """

    year = date.year
    month = date.month
    day = date.day
    h = date.hour
    m = date.minute
    s = date.second + date.microsecond * 1e-6

    jd = 367*year - int((7/4) * (year + int((month + 9)/12))) + int((275/9)*month) + day + 1721013.5
    ut = h/24 + m/(24*60) + s/(24*60**2)

    return jd + ut


def theta_g(jdn):
    """ Function to get Greenwich sidereal time from Julian Date
        This function follows Curtis eq. 5.50

        :param jdn: julian day number
        :return tg: Greenwich sidereal time (deg)
    """

    j0 = int(jdn) + 0.5
    ut = jdn - j0

    t0 = (j0 - j2000)/36525
    tg = 100.4606184 + 36000.77004*t0 + 0.000387933*t0**2 - 2.583e-8*t0**3

    tg += ut*360.98564724

    tg = functions.wrap_to_360(tg)

    return tg


def theta_g_2(jdn):
    # https://lweb.cfa.harvard.edu/~jzhao/times.html#ref7

    d = jdn - j2000
    T = d/36525

    gmst = 24110.54841 + 8640184.812866*T + 0.093104*T**2 - 0.0000062*T**3

    return gmst/60**2


def tle_to_datetime(epoch):
    """ Function to get a datetime object from a TLE epoch
        :param tle_epoch: epoch from tle
        :return date: datetime object of tle epoch
    """

    yr = int(epoch/1000)

    if yr < 57:
        year = yr + 2000
    else:
        year = yr + 1900

    dayfrac = epoch - yr*1000

    date = datetime.datetime(year, 1, 1) + datetime.timedelta(dayfrac - 1)

    return date


# def jd_to_date(jd):
#     """ Function to calculate Gregorian calendar date from Julian date
#
#         :param jd: julian date (days)
#         :return date: gregorian calendar date
#         :type date: tuple (day, month, year, hour, min, sec)
#     """
#
#     pass
#

# def jd_to_t0(jd, a=0):
#     """ Function to return Julian centuries since J2000 epoch
#         Follows Curtis Eqs. 5.49 - 5.52
#
#         :param jd: Julian date (days)
#         :param a: east longitude (deg)
#         :return T0: Julian centuries since J2000
#         :return theta_g: Greenwich sidereal time (deg)
#         :return theta: local sidereal time (deg)
#     """
#
#     jd2 = np.floor(jd) + 0.5
#     ut = jd - jd2
#     jd = jd2
#
#     T0 = (jd - j2000)/36525
#     theta_g = 100.4606184 + 36000.77004*T0 + 0.000387933*T0**2 - 2.583e-8*T0**3
#
#     theta_g = theta_g % 360 + 360.98564724*ut
#     theta = theta_g + a
#     theta = theta % 360
#
#     return T0, theta_g, theta
