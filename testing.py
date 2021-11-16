import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import dev_functions
import functions as func
import dev_functions as dev
from scipy.optimize import minimize

import maneuvers
import plotting
import propagators

def yline(val, max_x):
    plt.plot([0, max_x], [val, val], 'k--')


def xline(val, max_y):
    plt.plot([val, val], [0, max_y], 'k--')


def circle(r, cx, cy, style):
    ts = np.linspace(0, 2*np.pi, 100)
    xs = r*np.cos(ts) + cx
    ys = r*np.sin(ts) + cy

    plt.plot(xs, ys, style)


ra = 384648e3
soi = 66.1e6

mu_m = 4.9048695e12


def lunar_pc(r0, v0, phi0, lam1, rs=soi, rm=ra):
    # transfer orbit conditions
    E = v0**2/2 - cns.mu/r0
    h = v0*r0*np.cos(phi0)

    # patch point conditions
    r1 = np.sqrt(rs**2 + rm**2 - 2*rs*rm*np.cos(lam1))
    v1 = np.sqrt(2*(E + cns.mu/r1))
    phi1 = np.arccos(h/(r1*v1))
    gam1 = np.arcsin((rs/r1) * np.sin(lam1))

    # selenocentric orbit
    vm = np.sqrt(cns.mu/rm)
    v2 = np.sqrt(v1**2 + vm**2 - 2*v1*vm*np.cos(phi1 - gam1))
    eps2 = np.arcsin((vm/v2)*np.cos(lam1) - (v1/v2)*np.cos(lam1 + gam1 - phi1))

    h2 = v2*rs*abs(np.sin(eps2))
    E = v2**2/2 - mu_m/rs
    p = h2**2/mu_m

    e = np.sqrt(1 + 2*E*h2**2/mu_m**2)

    # perilune conditions
    rp = p/(1 + e)
    vp = h2/rp

    return rp, vp


# Initial orbit
r0 = cns.re + 200e3
v0 = np.sqrt(cns.mu/r0) + 3150
h0 = r0*v0
E0 = v0**2/2 - cns.mu/r0
a0 = -cns.mu/(2*E0)
e0 = np.sqrt(1 + 2*E0*h0**2/cns.mu**2)

# Patch point
lam1 = np.deg2rad(35)
r1 = np.sqrt(soi ** 2 + ra ** 2 - 2 * soi * ra * np.cos(lam1))
v1 = np.sqrt(2*(E0 + cns.mu/r1))

gam1 = np.arcsin((soi/r1) * np.sin(lam1))
theta = np.arccos((1/e0)*(h0**2/(cns.mu*r1) - 1))
phi1 = np.arccos(h0 / (r1 * v1))

# selenocentric orbit
vm = np.sqrt(cns.mu / ra)
v2 = np.sqrt(v1 ** 2 + vm ** 2 - 2 * v1 * vm * np.cos(phi1 - gam1))
eps2 = np.arcsin((vm / v2) * np.cos(lam1) - (v1 / v2) * np.cos(lam1 + gam1 - phi1))

h2 = v2 * soi * abs(np.sin(eps2))
E = v2 ** 2 / 2 - mu_m / soi
p = h2 ** 2 / mu_m
e = np.sqrt(1 + 2 * E * h2 ** 2 / mu_m ** 2)
a2 = -mu_m/(2*E)
e2 = e

theta2 = np.pi - np.arccos((1/e2)*(h2**2/(mu_m*soi) - 1))
asm = np.pi*2 - np.arccos(-1/e)

# perilune conditions
rp = p / (1 + e)
vp = h2 / rp

# orbit geometry
r, t, _, _ = func.perifocal_coords(a0, e0, np.linspace(0, np.pi*2, 10000))

t += gam1 + np.pi - theta
x = r*np.cos(t)
y = r*np.sin(t)
plt.plot(x, y)

# flyby geometry
r, t, _, _ = func.perifocal_coords(a2, e2)
print(360 - np.rad2deg(theta2))
print(360 - np.rad2deg(asm))

t += np.pi - theta2 - lam1
x = r*np.cos(t)
y = r*np.sin(t)
plt.plot(x - ra, y)

# patch point coords
xpp = -r1*np.cos(gam1)
ypp = -r1*np.sin(gam1)

# plot a bunch of geometry
plt.plot([0, -ra], [0, 0], 'k', linewidth=0.5)
plt.plot([0, xpp], [0, ypp], 'k', linewidth=0.5)
plt.plot([-ra, xpp], [0, ypp], 'k', linewidth=0.5)

# plot velocity vectors
v1_vec = np.array([-v1*np.sin(phi1 - gam1), -v1*np.cos(phi1 - gam1)])
vm_vec = np.array([0, vm])
v2_vec = v1_vec + vm_vec

plt.plot([xpp, xpp + vm_vec[0]*100e3], [ypp, ypp + vm_vec[1]*100e3])
plt.plot([xpp, xpp + v1_vec[0]*100e3], [ypp, ypp + v1_vec[1]*100e3])
plt.plot([xpp, xpp + v2_vec[0]*100e3], [ypp, ypp + v2_vec[1]*100e3])

# plot planets/SOIs
circle(cns.re, 0, 0, 'k')
circle(1737e3, -ra, 0, 'k')
circle(soi, -ra, 0, 'k--')
plt.axis('equal')
plt.show()
