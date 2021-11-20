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

    return np.cross(v1, v2)[2]


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
lam1 = np.radians(100)


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
v1_vec = np.array([-v1*np.sin(phi1 - gam1), -v1*np.cos(phi1 - gam1), 0])
r1_vec = np.array([-r1*np.cos(gam1), -r1*np.sin(gam1), 0])

a1, e1, _, _, w1, ta1 = func.vector_to_elements(r1_vec, v1_vec, b1.mu)


# Lunar flyby
rm1_vec = np.array([-b2.a, 0, 0])
vm_vec = np.array([0, -b1.v_circ(b2.a), 0])

r2_vec = r1_vec - rm1_vec
v2_vec = v1_vec - vm_vec
v2 = np.linalg.norm(v2_vec)
E2 = v2**2/2 - b2.mu/b2.soi
h2 = cross2d(r2_vec, v2_vec)

a2, e2, i2, _, w2, ta2 = func.vector_to_elements(r2_vec, v2_vec, b2.mu)

# w2 = 196
# Lunar escape
print(ta2)
ta3 = 360 - ta2
x3, y3, _ = func.elements_to_eci_pos(a2, e2, i2, 0, w2, ta3)
r3_vec = np.array([x3, y3, 0])
r4_vec = r3_vec + rm1_vec
vx, vy, _ = func.elements_to_eci_vel(a2, e2, 0, 0, w2, ta3, mu=b2.mu)
v3_vec = np.array([vx, vy, 0])
v4_vec = v3_vec + vm_vec


# plot shit
r, t, _, _ = func.perifocal_coords(a1, e1, np.linspace(0, 2*np.pi, 10000))
t += np.radians(w1)
x = r*np.cos(t)
y = r*np.sin(t)
plt.plot(x, y)

r, t, _, _ = func.perifocal_coords(a2, e2, np.linspace(np.radians(ta2 - 360), np.radians(360 - ta2), 1000))
t += np.radians(w2)
x = r*np.cos(t) + rm1_vec[0]
y = r*np.sin(t) + rm1_vec[1]
plt.plot(x, y)

scale = 25e3
vector(r1_vec[0], r1_vec[1], v1_vec, scale, 'bo-')
vector(r1_vec[0], r1_vec[1], v2_vec, scale, 'ro-')
vector(r1_vec[0], r1_vec[1], -vm_vec, scale, 'go-')

vector(r4_vec[0], r4_vec[1], v4_vec, scale, 'bo-')
vector(r4_vec[0], r4_vec[1], v3_vec, scale, 'ro-')
vector(r4_vec[0], r4_vec[1], vm_vec, scale, 'go-')

circle(b2.soi, rm1_vec[0], rm1_vec[1], 'k--')
circle(b1.r, 0, 0, 'k')
circle(b2.r, rm1_vec[0], rm1_vec[1], 'k')

plt.axis('equal')
plt.show()





















































# need to add a check to see if the selenocentric velocity is just going to result in leaving the SOI
# a, e, i, lan, w, ta = func.vector_to_elements(r1_vec, v1_vec, b1.mu)
#
# r, t, _, _ = func.perifocal_coords(a, e, np.linspace(np.radians(ta - 360), np.radians(360 - ta), 1000))
# t += np.radians(w)
#
#
# x = r*np.cos(t)
# y = r*np.sin(t)
# plt.plot(x, y)
#
# a, e, i, lan, w, ta = func.vector_to_elements(r2_vec, v2_vec, b2.mu)
#
# # w = 180
# print(w)
# if h2 < 0:
#     w = 360 - w
#
# r, t, _, _ = func.perifocal_coords(a, e, np.linspace(np.radians(ta - 360), np.radians(360 - ta), 1000))
# t += np.radians(w)
#
#
# x = r*np.cos(t) + rm1_vec[0]
# y = r*np.sin(t) + rm1_vec[1]
#
# circle(b2.soi, rm1_vec[0], rm1_vec[1], 'k--')
#
# plt.plot(x, y)
#
# vector(0, 0, rm1_vec, style='ko-')
# vector(0, 0, r1_vec, style='ko-')
#
# vector(r1_vec[0], r1_vec[1], v1_vec, 25e3, 'bo-')
# vector(r1_vec[0], r1_vec[1], v2_vec, 25e3, 'ro-')
# vector(r1_vec[0], r1_vec[1], vm_vec, 25e3, 'go-')
#
# plt.axis('equal')
# plt.show()
#
#
#
# if h2 > 0:
#     # Prograde orbit
#     pass
# elif h2 < 0:
#     # retrograde orbit
#     pass
# else:
#     # radial orbit
#     pass
#
# # print(v1)
# # print(np.linalg.norm(v2_vec))
# #
# # vector(0, 0, v1_vec, style='bo-')
# # vector(0, 0, vm_vec, style='go-')
# # vector(0, 0, v2_vec, style='ro-')
# # vector(0, 0, r2_vec)
# # plt.show()
