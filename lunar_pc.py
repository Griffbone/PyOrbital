import numpy as np
import constants as cns
import functions as func
import matplotlib.pyplot as plt


def circle(r, cx, cy, style):
    ts = np.linspace(0, 2*np.pi, 100)
    xs = r*np.cos(ts) + cx
    ys = r*np.sin(ts) + cy

    plt.plot(xs, ys, style)


class body:
    def __init__(self, r, mu, soi, a):
        """ Initialize body class
            :param r: body radius
            :param mu: body gravitational parameter
            :param a: orbit semimajor axis
        """

        self.r = r
        self.mu = mu
        self.soi = soi
        self.a = a

    def v_circ(self, r):
        """ Get circular orbit velocity at given radius
            :param r: orbit radius
        """

        return np.sqrt(self.mu/r)

    def v_esc(self, r):
        """ Get escape velocity at given radius
            :param  r: orbit radius
        """

        return np.sqrt(2*self.mu/r)


earth = body(6378e3, cns.mu, 0, 0)
moon = body(1737447.78, 4.9048695e12, 66.1e6, 384648e3)

b1 = earth
b2 = moon

# Inputs (2 initial, one patch point)
r0 = b1.r + 200e3
v0 = b1.v_circ(r0) + 3150
lam1 = np.radians(37.25)

# Initial orbit
h0 = r0*v0
E0 = v0**2/2 - b1.mu/r0
a0 = -b1.mu/(2*E0)
e0 = np.sqrt(1 + 2*E0*h0**2/b1.mu**2)

# Patch point
r1 = np.sqrt(b2.soi**2 + b2.a**2 - 2*b2.soi*b2.a*np.cos(lam1))
v1 = np.sqrt(2*(E0 + b1.mu/r1))
gam1 = np.arcsin((b2.soi/r1)*np.sin(lam1))
ta1 = np.arccos((1/e0)*(h0**2/(b1.mu*r1) - 1))
phi1 = np.arccos(h0/(r1*v1))

# Lunar orbit
vm = np.sqrt(b1.mu/b2.a)

v2 = np.sqrt(v1**2 + vm**2 - 2*v1*vm*np.cos(phi1 - gam1))
eps = np.arcsin((vm/v2)*np.cos(lam1) - (v1/v2)*np.cos(lam1 + gam1 - phi1))

h2 = v2*b2.soi*abs(np.sin(eps))
E2 = v2**2/2 - b2.mu/b2.soi
p2 = h2**2/b2.mu
e2 = np.sqrt(1 + 2*E2*h2**2/b2.mu**2)
a2 = -b2.mu/(2*E2)
ta2 = 2*np.pi - np.arccos((1/e2)*(h2**2/(b2.mu*b2.soi) - 1))

if eps > 0:
    w2 = -ta2 - lam1
else:
    w2 = ta2 - lam1

# Escape patch point
nu1 = ta2
nu2 = 2*np.pi - nu1

F = np.arctanh(np.sqrt((e2 - 1)/(e2 + 1))*np.tan(nu2/2))*2
M = e2*np.sinh(F) - F
tof = ((h2**3/b2.mu**2)/(e2**2 - 1)**(3/2))*M*2
dnu = 2*np.pi*(tof/func.period(b2.a, b1.mu))


