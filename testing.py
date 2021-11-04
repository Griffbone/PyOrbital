import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import functions as func
import dev_functions as dev
from scipy.optimize import minimize


def circle(r, cx, cy, style):
    ts = np.linspace(0, 2*np.pi, 100)
    xs = r*np.cos(ts) + cx
    ys = r*np.sin(ts) + cy

    plt.plot(xs, ys, style)


def vector_angle(v1, v2):
    dot = np.dot(v1, v2)
    mag = np.linalg.norm(v1)*np.linalg.norm(v2)

    return np.arccos(dot/mag)


def patch_point(sat, targ, rs):
    if targ.mu != sat.mu:
        return None
    elif sat.ap < targ.pe - rs:
        return None
    elif sat.pe > targ.ap + rs:
        return None

    ta1 = np.degrees(np.arccos(sat.a*(1 - sat.e**2)/((targ.a - rs)*sat.e) - 1/sat.e))
    ta2 = 360 - ta1

    def dist(ta_sat):
        Dt = func.t_between(sat.a, sat.e, sat.ta, ta_sat, sat.mu)
        ta_targ = func.ta_change(targ.a, targ.e, targ.ta, Dt, targ.mu)

        print(ta_sat)
        x, y, z = sat.get_eci(ta_sat)
        r1 = np.array([x, y, z])

        x, y, z = targ.get_eci(ta_targ)
        r2 = np.array([x, y, z])

        D = np.linalg.norm(r1 - r2)

        return abs(D - rs)

    ta = minimize(dist, (ta1 - ta2)/2, bounds=(ta1, ta2))


def dist(sat, targ, ta_sat, rs):
    Dt = func.t_between(sat.a, sat.e, sat.ta, ta_sat, sat.mu)
    ta_targ = func.ta_change(targ.a, targ.e, targ.ta, Dt, targ.mu)

    # print(ta_sat)
    x, y, z = sat.get_eci(ta_sat)
    r1 = np.array([x, y, z])

    x, y, z = targ.get_eci(ta_targ)
    r2 = np.array([x, y, z])

    D = np.linalg.norm(r1 - r2)

    return abs(D - rs), Dt


ra = 384000e3
rp = 200000 + cns.re
rs = 66.1e6

sat = func.Elements((ra + rp)/2, (ra - rp)/(ra + rp), 0, 0, 0, 0, cns.mu)
moon = func.Elements(0.3844e9, 0, 0, 0, 0, 150, cns.mu)

#############################################################################
ts = np.linspace(0, 360*2, 1000)
Ds = []
dts = []
times = []

time = 0

print(sat.a)
print(sat.e)

for t in ts:
    time = func.t_between(sat.a, sat.e, sat.ta, t, sat.mu)
    times.append(time)

# plt.plot(ts, Ds)
# plt.show()
plt.plot(ts, times)
plt.show()
