from pyOrbital import constants as cons
from numpy import sin, cos, sqrt


def thrust(t):
    return 5e3

def drag(y, v):
    return (1/2)*v^2

def mass(time):
    return 1e3

def gravity(y):
    return cons.mu/(y + cons.re)**2

def gravity_turn_odes(t, y):
    x, y, v, gam = y

    xdot = v*cos(gam)
    ydot = v*sin(gam)
    vdot = thrust(t)/mass(t) - drag(y, v) - gravity(y)*sin(gam)
    gamdot = -gravity(y)*cos(gam)/v + ((v**2)/(cons.re + y))*cos(gam)

    #
    # v, h, gam =

print(gravity(0.0))