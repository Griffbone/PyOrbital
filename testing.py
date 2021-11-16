import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import dev_functions
import functions as func
import dev_functions as dev
from scipy.optimize import minimize
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


ra = 0.3744e9
rp = 200000 + cns.re
soi = 66.1e6

sat = func.Elements((ra + rp)/2, (ra - rp)/(ra + rp), 0, 0, 0, 0, cns.mu)
moon = func.Elements(0.3844e9, 0, 0, 0, 0, 140, cns.mu)

mu_m = 4.9028e3*1000**3


def lunar_pc(r0, v0, phi0, lam1, rs=soi, rm=ra):
    # transfer orbit conditions
    E = v0**2/cns.mu - cns.mu/r0**2
    h = v0*r0*np.sin(phi0)

    # patch point conditions
    r1 = np.sqrt(rs**2 + rm**2 - 2*rs*rm*np.cos(lam1))
    v1 = np.sqrt(2*(E + cns.mu/r1))
    phi1 = np.arccos(h/(r1*v1))
    gam1 = np.arcsin((rs/r1) * np.sin(lam1))

    # selenocentric orbit
    vm = np.sqrt(cns.mu/rm)
    v2 = np.sqrt(v1**2 + vm**2 - 2*v1*vm*np.cos(phi1 - gam1))
    eps2 = np.sin((vm/v2)*np.cos(lam1) - (v1/v2)*np.cos(lam1 + gam1 - phi1))

    h2 = v2*rs*np.sin(eps2)
    E = v2**2/2 - mu_m/rs
    p = h2**2/mu_m

    e = np.sqrt(1 + 2*E*h**2/mu_m**2)
    rp = p/(1 + e)

    print(rp)


print(soi)
print(ra)

lunar_pc(6697e3, 10.8462e3, 0, np.deg2rad(30))
