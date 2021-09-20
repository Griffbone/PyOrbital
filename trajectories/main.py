import scipy as sp
import numpy as np
import scipy.integrate
import scipy.optimize as opt
import matplotlib.pyplot as plt

global g, re, T
g = 9.80665
re = 6371e3


def ydot(t, y, a, b):
    x, y, vx, vy, m = y

    gam = np.arctan2(vy, vx)

    theta = np.arctan(a*t + b)

    if vx**2/(re + y) - g == 0:
        T = 0
    else:
        T = 100000 * g * 0.6

    ax = T*np.cos(theta)/m
    ay = (T*np.sin(theta) - g*m + m*vx**2/(re + y))/m

    return [vx, vy, ax, ay, -T/(400*g)]


m0 = 100000
# T = m0*g*1.2
# tmax = 90000/(T/(400*g))
tmax = 300
y0 = [0, 100e3, 400, 400, m0]


    # def h_fun(vars):
    #     a, b = vars
    #
    #     solution = scipy.integrate.solve_ivp(lambda t, y: ydot(t, y, a, b), [0, tmax], y0, dense_output=True, t_eval=np.linspace(0, tmax, 1000))
    #
    #     g1 = 1/1000
    #     g2 = 100
    #
    #     gam_err = abs(np.arctan2(solution.y[3][-1], solution.y[2][-1]))     * g2
    #     a_err = abs(200e3 - (solution.y[1][-1]))                            * g1
    #     H = a_err + gam_err
    #
    #     return H
    #
    # sol = opt.minimize(h_fun, [-0.1, 10])
    #
    # a, b = sol.x
    #
    # print([a, b])

a = 0
b = 1000
solution = scipy.integrate.solve_ivp(lambda t, y: ydot(t, y, a, b), [0, tmax], y0, dense_output=True, t_eval=np.linspace(0, tmax, 1000))

t = solution.t
y = solution.y
xs = y[0]
ys = y[1]
vx = y[2]
vy = y[3]
ms = y[4]

plt.subplot(1, 2, 1)
plt.plot(xs, ys)
plt.axis('equal')

plt.subplot(1, 2, 2)
plt.plot(t, vx)
plt.plot(t, vy)
plt.legend(['vx', 'vy'])

# plt.figure()
# rho = 1.225*np.exp(-g*0.02896648*ys/(8.3144*300))
# print(max(rho))
# # print(np.shape(rho))
#
# q = [np.sqrt(v**2 + u**2) for v, u in zip(vx, vy)]
# q = 0.5*rho*np.array(q)**2
# plt.plot(t, q/101325)

plt.show()
