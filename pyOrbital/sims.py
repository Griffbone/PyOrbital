import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import RK45, solve_ivp, solve_bvp
import pyOrbital.functions as func
import pyOrbital.constants as const


def plot_planet(r):
    ts = np.linspace(0, 2*np.pi, 1000)
    x = r*np.cos(ts)
    y = r*np.sin(ts)
    plt.fill(x, y, 'b', alpha=0.5, edgecolor='k')


def acc_body(ro):
    """ Function to calculate gravity vector towards large central body (M >> m)
        :param ro: position of the child body
    """

    # Acceleration due to gravity
    mu = const.mu
    roc = -ro
    r = np.linalg.norm(roc)
    a = mu/r**2

    # J2 perturbation acceleration
    J2 = 1.082629e-3
    x = ro[0]
    y = ro[1]
    z = ro[2]

    v = np.array([
        (1 - (3*z**2)/r**2)*(x/r),
        (1 - (3*z**2)/r**2)*(y/r),
        (3 - (3*z**2)/r**2)*(z/r)
    ])

    c = (-3/2)*((mu*J2*const.re**2)/r**4)

    # Return complete acceleration
    return a*(roc/r)# + c*v


def derivatives(t, y):
    ro = y[0:3]
    vc = y[3:]

    rdot = vc
    vdot = acc_body(ro)

    return np.concatenate([rdot, vdot])


J2 = 1.082629e-3
mu = const.mu

r = const.re + 397000
incl = np.radians(0)
n = 5

T = 2*np.pi*np.sqrt(r**3/mu)
r = r*np.array([1, 0, 0])
v = np.sqrt(mu/np.linalg.norm(r))
v = v*np.array([0, np.cos(incl), np.sin(incl)])
# a, e, i, omega, w, M = func.vector_to_elements(r, v, mu)
#
# dodt = -(3/2)*J2*(6371e3/(a*(1 - e**2)))**2*(mu/a**3)*np.cos(i)

y0 = np.concatenate([r, v]) #DOP853
out = solve_ivp(derivatives, [0, T*n], y0, t_eval=np.linspace(0, T*n, 1000000), method='DOP853')
x, y, z, vx, vy, vz = out.y

saxs = []
es = []
incls = []
omegas = []

# for i in range(0, len(x)):
#     r = np.array([x[i], y[i], z[i]])
#     v = np.array([vx[i], vy[i], vz[i]])
#
#     a, e, i, omega, w, M = func.vector_to_elements(r, v, mu)
#
#     saxs.append(a)
#     es.append(e)
#     incls.append(i)
#     omegas.append(omega)


# plt.plot(out.t/3600/24, np.degrees(incls))
# plt.plot(out.t/3600/24, np.degrees(omegas))
# plt.plot(out.t/3600/24, np.array(saxs)/max(saxs))
# plt.legend(['Inclination', 'RAAN'])
# plt.show()

# ax = plt.axes(projection='3d')
plot_planet(const.re/1000)
plt.plot(x/1000, y/1000)#, z/1000)
plt.axis('equal')
plt.show()
