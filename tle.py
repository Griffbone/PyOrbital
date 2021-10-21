import re
import constants as cns
import numpy as np
import astrotime as at


def read_l2(line):
    words = re.split('\s', line)
    words = [w for w in words if w != '']

    num_1 = words[1]
    desig = words[2]
    epoch = words[3]
    ndot = float(words[4])

    epoch = re.split('\.', epoch)
    frac = epoch[1]
    day = epoch[0][2::]
    year = epoch[0][0:2]

    jd = cns.j2000 + 365.25*float(year) + float(day) + float(frac)*10**(-len(frac))

    nddot = re.split('-', words[5])
    nddot = float(nddot[0]) * 10**(-float(nddot[1]))

    drag = re.split('-', words[6])
    drag = float(drag[0]) * 10**(-float(drag[1]))

    typ = int(words[7])
    csum = int(words[8])

    return num_1, desig, epoch, ndot, nddot, drag, typ, csum, jd


def read_l3(line):
    words = re.split('\s', line)
    words = [w for w in words if w != '']

    num_2 = words[1]
    inc = float(words[2])
    raan = float(words[3])
    ecc = float(words[4]) * 10**(-7)
    w = float(words[5])
    ma = float(words[6])
    n = float(words[7])

    return num_2, inc, raan, ecc, w, ma, n


def read_tle(l1, l2, l3):
    name = l1
    num_1, desig, epoch, ndot, nddot, drag, typ, csum, jd = read_l2(l2)
    num_2, inc, raan, ecc, w, ma, n = read_l3(l3)

    n_rad = n*(1/(60**2*24))*2*np.pi
    sma = (n_rad**2/cns.mu)**(-1/3)

    tle = {
        'name': name,
        'number': num_1,
        'designation': desig,
        'epoch': epoch,
        'ndot': ndot,
        'nddot': nddot,
        'drag': drag,
        'ephm type': typ,
        'checksum': csum,
        'inclination': inc,
        'raan': raan,
        'eccentricity': ecc,
        'argpe': w,
        'ma': ma,
        'n': n,
        'sma': sma,
        'jd': jd
    }

    return tle


def epoch_to_jd(epoch):
    year = int(epoch[0:2])

    if year < 57:
        year += 2000
    else:
        year += 1900
