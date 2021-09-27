import scipy as sp
import numpy as np
import scipy.integrate
import scipy.optimize as opt
import matplotlib.pyplot as plt
import pyOrbital.functions as fncs
import pyOrbital.constants as cons

global g, re, T
g = 9.80665
re = 6371e3


def flat_earth_to_ap_pe(y, vx, vy):
    r = np.array([cons.re + y, 0, 0])
    v = np.array([vy, vx, 0])

    h = np.cross(r, v)
    evec = np.cross(v, h)/cons.mu - r/np.linalg.norm(r)
    e = np.linalg.norm(evec)

    a = 1/((2/np.linalg.norm(r)) - (np.linalg.norm(v)**2/cons.mu))

    rp = a*(1 - e)
    ra = a*(1 + e)

    return rp - cons.re, ra - cons.re


def ydot(t, y, a, b):
    x, y, vx, vy, m = y

    theta = np.arctan(a*t + b)

    if t <= 800:
        T = 100e3 * g * 0.6
    else:
        T = 0

    ax = T*np.cos(theta)/m
    ay = (T*np.sin(theta) - g*m + m*vx**2/(re + y))/m

    return [vx, vy, ax, ay, -T/(400*g)]


m0 = 100000
y0 = [0, 100e3, np.sqrt(2)*1500, np.sqrt(2)*1500, m0]


def h_fun(vars):
    tmax, a, b = vars

    solution = scipy.integrate.solve_ivp(lambda t, y: ydot(t, y, a, b), [0, tmax], y0, dense_output=True, t_eval=np.linspace(0, tmax, 1000))

    y = solution.y
    xf = y[0][-1]
    yf = y[1][-1]

    vxf = y[2][-1]
    vyf = y[3][-1]

    pe, ap = flat_earth_to_ap_pe(yf, vxf, vyf)

    return (pe - 200e3)**2 + (ap - 200e3)**2


sol = opt.minimize(h_fun, np.array([250, -0.1, 10]))

tmax, a, b = sol.x
solution = scipy.integrate.solve_ivp(lambda t, y: ydot(t, y, a, b), [0, tmax], y0, dense_output=True, t_eval=np.linspace(0, tmax, 1000))

t = solution.t
y = solution.y
xs = y[0]
ys = y[1]
vx = y[2]
vy = y[3]
ms = y[4]
thetas = [np.arctan(a*x + b) for x in t]

plt.subplot(1, 2, 1)
plt.plot(xs/1e3, ys/1e3)
plt.axis('equal')

plt.subplot(1, 2, 2)
plt.plot(t, vx)
plt.plot(t, vy)
plt.legend(['vx', 'vy'])

plt.show()

vg = scipy.integrate.trapz(g*np.sin(thetas), t)
print(vg + np.sqrt(cons.mu/(cons.re + 200e3)))

# plt.plot(t, ms)
# plt.show()
# plt.figure()
# rho = 1.225*np.exp(-g*0.02896648*ys/(8.3144*300))
# print(max(rho))
# # print(np.shape(rho))
#
# q = [np.sqrt(v**2 + u**2) for v, u in zip(vx, vy)]
# q = 0.5*rho*np.array(q)**2
# plt.plot(t, q/101325)
