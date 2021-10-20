import re

l1 = 'ISS (ZARYA)'
l2 = '1 25544U 98067A   21292.89044074  .00004617  00000-0  92733-4 0  9998'
l3 = '2 25544  51.6431  86.1166 0004166 125.7553 304.8177 15.48743137307929'


def read_l1(line):
    words = re.split('\s', line)
    words = [w for w in words if w != '']

    id = words[0]
    line = words[1][words[1].find('(')+1:words[1].find(')')]

    return id, line


def read_l2(line):
    words = re.split('\s', line)
    words = [w for w in words if w != '']

    num_1 = words[0]
    desig = words[2]

    print(words)

read_l2(l2)
