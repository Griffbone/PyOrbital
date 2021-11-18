import numpy as np
import constants as cns
import functions as func
import matplotlib.pyplot as plt


def circle(r, cx, cy, style):
    ts = np.linspace(0, 2*np.pi, 100)
    xs = r*np.cos(ts) + cx
    ys = r*np.sin(ts) + cy

    plt.plot(xs, ys, style)


def vector(cx, cy, vec, scale=1, style='ko-'):
    plt.plot([cx, cx + vec[0]*scale], [cy, cy + vec[1]*scale], style, ms=3, markevery=[-1])


def cross2d(v1, v2):
    v1 = np.array([v1[0], v1[1], 0])
    v2 = np.array([v2[0], v2[1], 0])

    return np.linalg.norm(np.cross(v1, v2))


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
phi1 = np.arccos(h0/(r1*v1))
v1_vec = np.array([-v1*np.sin(phi1 - gam1), -v1*np.cos(phi1 - gam1)])

# Lunar flyby
vm_vec = np.array([0, -b1.v_circ(b2.a)])
v2_vec = v1_vec - vm_vec
v2 = np.linalg.norm(v2_vec)
E2 = v2**2/2 - b2.mu/b2.soi
h2 = cross2d(v2, )

print(v1)
print(np.linalg.norm(v2_vec))

vector(0, 0, v1_vec, style='bo-')
vector(0, 0, vm_vec, style='go-')
vector(0, 0, v2_vec, style='ro-')
plt.show()