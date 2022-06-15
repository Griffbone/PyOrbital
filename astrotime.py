"""
    Functions for conversion between time systems.
    Author      : Griffin Jourda
    Date        : 2/16/22
    Last Edited : 2/16/22

    Functions
        date_to_jd  :   convert gregorian date to Julian date
        jd_to_gmst  :   calculate Greenwich mean sidereal time from UT1 JD
        jd_to_lst   :   calculate local sidereal time from UT1 JD and longitude
        deg_to_hms  :   convert angle in degrees to hours:minutes:seconds format
        hms_to_deg  :   convert hours:minutes:seconds format to angle in degrees

    Dependencies
        constants       :   j2000
        numpy           :   pi, deg2rad
        math_functions  :   wrap_to_2pi
"""


import constants as cns
import numpy as np
import math_functions as mfun


def date_to_jd(year, month, day, h, m, s):
    """ Convert Gregorian calendar date to Julian day number and Julian day fraction.
        Follows Curtis Eq. 5.48. Valid for March 1 1900 to Feb 28 2100.

        :param year: gregorian calendar year
        :param month: gregorian calendar month
        :param day: gregorian calendar day
        :param h: hours
        :param m: minutes
        :param s: seconds

        :return jdn: Julian day number (days since Julian calendar epoch)
        :return jd_frac: Julian day fraction (days since jdn)
    """

    jdn = 367*year - int((7/4) * (year + int((month + 9)/12))) + int((275/9)*month) + day + 1721013.5
    jd_frac = s/86400 + m/1440 + h/24

    if jd_frac > 1:
        jdn += int(jd_frac)
        jd_frac = jd_frac - int(jd_frac)

    return jdn, jd_frac


def jd_to_day(jd):
    """ Function to convert Julian date to Gregorian Calendar date
        :param jd: Julian date

        :return y: year
        :return m: month
        :return d: day
        :return h: hour
        :return m: minute
        :return s: second
    """
    T1900 = (jd - 2415019.5)/365.25
    year = 1900 + int(T1900)
    leapyrs = int((year - 1900 - 1)*0.25)
    days = (jd - 2415019.5) - ((year - 1900)*365 + leapyrs)

    lmonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    # Check if we are actually in the previous year
    if days < 1:
        year = year - 1
        leapyrs = int((year - 1900 - 1)*0.25)
        days = (jd - 2415019.5) - ((year - 1900)*365 + leapyrs)

    # Check if this year is a leap year
    if (year % 4) == 0:
        lmonth[1] = 29

    dayofyr = int(days)
    ndays = 0
    month = 0

    while ndays < dayofyr:
        ndays += lmonth[month]
        month += 1

    day = dayofyr - sum(lmonth[0:month - 1])
    tau = 24*(days - dayofyr)
    h = int(tau)
    m = int((tau - h)*60)
    s = (tau - h - m/60)*3600

    return year, month, day, h, m, s


def j2000_cents(jd):
    """ Calculate centuries elapsed since J2000 epoch from Julian date

        :param jd: Julian date
        :return T: Centuries since J2000 epoch
    """

    T = (jd - cns.j2000)/36525
    return T


def jd_to_gmst(jd_ut1):
    """ Calculate Greenwich mean sidereal time from UT1 Julian date.
        Function follows Vallado Algorithm 15. Valid for years 1901 to 2199.

        :param jd_ut1: Julian day in UT1 time
        :return t_gmst: Greenwich mean sidereal time (rad)
    """

    T = (jd_ut1 - cns.j2000)/36525
    t_gmst = 67310.54841 + (876600*60*60 + 8640184.812866)*T + 0.093104*T**2 - 6.2e-6*T**3
    t_gmst *= (1/240)*(np.pi/180)
    t_gmst = mfun.wrap_to_2pi(t_gmst)

    return t_gmst


def jd_to_lst(jd_ut1, lon):
    """ Calculate local sidereal time from UT1 Julian date and longitude.
        Valid for years 1901 to 2199.

        :param jd_ut1: Julian day in UT1 time
        :param lon: east longitude (deg); west longitudes entered as negative
        :return lst: local sidereal time (rad)
    """

    lst = jd_to_gmst(jd_ut1) + np.deg2rad(lon)
    lst = mfun.wrap_to_2pi(lst)

    return lst


def deg_to_hms(theta):
    """ Convert an angle from degrees to hours:minutes:seconds

        :param theta: angle to convert (deg)
        :return h: hours
        :return m: minutes
        :return s: seconds
    """

    h = int(theta/15)
    theta -= h*15

    m = int(theta/0.25)
    theta -= m*0.25

    s = theta*240

    return h, m, s


def hms_to_deg(h, m, s):
    """ Convert an angle from hours:minutes:seconds to degrees.

        :param h: hours
        :param m: minutes
        :param s: seconds
        :return theta: angle (deg)
    """

    theta = h*15 + m*0.25 + s/240
    return theta
