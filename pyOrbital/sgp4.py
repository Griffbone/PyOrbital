import numpy as np
import re


def parse_tle(fname):
    lines = []
    with open(fname) as f:
        n = 0
        for line in f:
            line = line.split(' ')

            newline = []
            for i in line:
                i = i.split('\n')[0]
                if i != '' and i != '\n':
                    newline.append(i)

            lines.append(newline)

    l1 = lines[1]
    l2 = lines[2]

    epoch = l1[3]
    ndot = l1[4]
    nddot = l1[5]
    drag = l1[6]
    ephm = l1[7]
    csum = l1[8]

    incl = l2[2]
    raan = l2[3]
    ecc = l2[4]
    argpe = l2[5]
    ma = l2[6]
    n = l2[7]
    rev = l2[8]

    print([incl, raan, ecc, argpe])



parse_tle('ISS.txt')


