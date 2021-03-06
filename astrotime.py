def date_to_jd(day, month, year, h, m, s):
    """ Function to calculate Julian day number from Gregorian Calendar date
        Follows Curtis Eq. 5.48

        :param day: gregorian calendar day
        :param month: gregorian calendar month
        :param year: gregorian calendar year
        :return jdn: Julian day number (days)
    """
    # Valid range: March 1 1900 - Feb 28 2100

    jd = 367*year - int((7/4) * (year + int((month + 9)/12))) + int((275/9)*month) + day + 1721013.5
    ut = ((s/60 + m)/60 + h)/24

    jdn = jd + ut

    return jdn


def gmst(jdn):
    """ Function to calculate Greenwich mean sidereal time from Julian day number
        Function follows Vallado Algorithm 15

        :param jdn: julian day number
        :return t_gmst: Greenwich mean sidreal time (deg)
    """
    # Valid range: 1901 - 2199
    T = (jdn - 2451545.0) / 36525
    t_gmst = 67310.54841 + (876600 * 60 * 60 + 8640184.812866) * T + 0.093104 * T ** 2 - 6.2e-6 * T ** 3

    t_gmst -= int(t_gmst / 86400) * 86400
    t_gmst *= 1 / 240

    if t_gmst < 0:
        t_gmst += 360

    return t_gmst


def lst(jdn, lon):
    """ Function to calculate local sidereal time from julian day number
        :param jdn: Julian day number
        :param lon: east longitude (deg)
        :return lst: local sidereal time (deg)
    """

    lst = gmst(jdn) + lon
    lst -= int(lst/360)*360

    if lst < 0:
        lst += 360

    return lst
