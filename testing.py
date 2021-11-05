import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import dev_functions
import functions as func
import dev_functions as dev
from scipy.optimize import minimize


import astrotime as at

r = np.array([-6045, -3490, 2500])*1000
v = np.array([-3.457, 6.618, 2.533])*1000

print(func.vector_to_elements(r, v, cns.mu))










# sat.perifocal_plot()

# plt.show()


# import propagators
#
#
# def yline(val, max_x):
#     plt.plot([0, max_x], [val, val], 'k--')
#
#
# def xline(val, max_y):
#     plt.plot([val, val], [0, max_y], 'k--')
#
#
# def circle(r, cx, cy, style):
#     ts = np.linspace(0, 2*np.pi, 100)
#     xs = r*np.cos(ts) + cx
#     ys = r*np.sin(ts) + cy
#
#     plt.plot(xs, ys, style)
#
#
# def vector_angle(v1, v2):
#     dot = np.dot(v1, v2)
#     mag = np.linalg.norm(v1)*np.linalg.norm(v2)
#
#     return np.arccos(dot/mag)
#
#
# def patch_point(sat, targ, rs):
#     if targ.mu != sat.mu:
#         return None
#     elif sat.ap < targ.pe - rs:
#         return None
#     elif sat.pe > targ.ap + rs:
#         return None
#
#     ta1 = np.degrees(np.arccos(sat.a*(1 - sat.e**2)/((targ.a - rs)*sat.e) - 1/sat.e))
#     ta2 = 360 - ta1
#
#     def dist(ta_sat):
#         Dt = func.t_between(sat.a, sat.e, sat.ta, ta_sat, sat.mu)
#         ta_targ = func.ta_change(targ.a, targ.e, targ.ta, Dt, targ.mu)
#
#         print(ta_sat)
#         x, y, z = sat.get_eci(ta_sat)
#         r1 = np.array([x, y, z])
#
#         x, y, z = targ.get_eci(ta_targ)
#         r2 = np.array([x, y, z])
#
#         D = np.linalg.norm(r1 - r2)
#
#         return abs(D - rs)
#
#     ta = minimize(dist, (ta1 - ta2)/2, bounds=(ta1, ta2))
#
#
# def dist(sat, targ, ta_sat, rs):
#     Dt = func.t_between(sat.a, sat.e, sat.ta, ta_sat, sat.mu)
#     ta_targ = func.ta_change(targ.a, targ.e, targ.ta, Dt, targ.mu)
#
#     # print(ta_sat)
#     x, y, z = sat.get_eci(ta_sat)
#     r1 = np.array([x, y, z])
#
#     x, y, z = targ.get_eci(ta_targ)
#     r2 = np.array([x, y, z])
#
#     D = np.linalg.norm(r1 - r2)
#
#     return abs(D - rs), Dt
#
#
# ra = 0.3744e9
# rp = 200000 + cns.re
# soi = 66.1e6
#
# sat = func.Elements((ra + rp)/2, (ra - rp)/(ra + rp), 0, 0, 0, 0, cns.mu)
# moon = func.Elements(0.3844e9, 0, 0, 0, 0, 140, cns.mu)
#
# #############################################################################
# ts, xm, ym, zm = propagators.kepler_propagation(moon.a, moon.e, moon.i, 0, 0, moon.ta, sat.T/2)
# _, xs, ys, zs = propagators.kepler_propagation(sat.a, sat.e, sat.i, 0, 0, sat.ta, max(ts))
# rs = np.sqrt((xs - xm)**2 + (ys - ym)**2 + (zs - zm)**2)
#
# tas, tam, _ = (dev_functions.patch_point(sat, moon, 66.1e6))
# Dt = func.t_between(sat.a, sat.e, sat.ta, tas, sat.mu)
#
# print(tas)
#
# x, y, _ = sat.get_eci(tas)
# plt.scatter(x, y)
#
# x, y, _ = moon.get_eci(tam)
# plt.scatter(x, y)
# circle(soi, x, y, 'k--')
#
# sat.perifocal_plot()
# moon.perifocal_plot()
#
# plt.show()
