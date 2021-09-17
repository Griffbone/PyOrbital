import numpy as np
import functions as func

days = 60*60*24


def stumpff_s(z):
    if z > 0:
        s = (np.sqrt(z) - np.sin(np.sqrt(z)))/(np.sqrt(z))**3
    elif z < 0:
        s = (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/(np.sqrt(-z))**3
    else:
        s = 1/6

    return s


def stumpff_c(z):
    if z > 0:
        c = (1 - np.cos(np.sqrt(z)))/z
    elif z < 0:
        c = (np.cosh(np.sqrt(-z)) - 1)/(-z)
    else:
        c = 0.5

    return c


def lambert(r1, r2, mu, dt):
    """ Solve Lambert's problem with method laid out in Curtis section 5.3
        :param r1: position vector 1
        :param r2: position vector 2
        :param dt: time of flight between positions one and two
        :return v1: velocity at initial position
        :return v2: velocity at final position
    """

    r1_m = np.linalg.norm(r1)
    r2_m = np.linalg.norm(r2)

    # Determine true anomaly change
    dta = np.arccos(np.dot(r1, r2)/(r1_m*r2_m))
    z = np.cross(r1, r2)[2]
    cosi = z/(r1_m*r2_m*np.sin(dta))

    if cosi >= 0:
        if z >= 0:
            dta = np.arccos(np.dot(r1, r2)/(r1_m*r2_m))
        else:
            dta = np.pi*2 - np.arccos(np.dot(r1, r2)/(r1_m*r2_m))
    elif cosi < 0:
        if z < 0:
            dta = np.arccos(np.dot(r1, r2)/(r1_m*r2_m))
        else:
            dta = np.pi*2 - np.arccos(np.dot(r1, r2)/(r1_m*r2_m))

    # Apply method of universal variables
    A = np.sin(dta)*np.sqrt((r1_m*r2_m)/(1 - np.cos(dta)))

    def F(z):
        S = stumpff_s(z)
        C = stumpff_c(z)

        y = r1_m + r2_m + A*((z*S - 1)/(np.sqrt(C)))

        return ((y/C)**(3/2))*S + A*np.sqrt(y) - np.sqrt(mu)*dt

    z = fsolve(F, np.array([1.5]))
    z = z[0]

    y = r1_m + r2_m + A * ((z * stumpff_s(z) - 1) / (np.sqrt(stumpff_c(z))))

    # Lagrange coefficients
    f = 1 - y/r1_m
    g = A*np.sqrt(y/mu)
    gdot = 1 - y/r2_m

    # Initial and final velocities
    v1 = (1/g)*(r2 - f*r1)
    v2 = (1/g)*(gdot*r2 - r1)

    return v1, v2


def gibbs(r1, r2, r3, mu, tol=10e-3):
    """ Function to solve Gibb's problem of preliminary orbit determination
        Based on the algorithm laid on in Curtis section 5.2

        :param r1: first position vector (m)
        :param r2: second position vector (m)
        :param r3: third position vector (m)
        :param mu: gravitational parameter
        :param tol: co-planar tolerance (dot product quantity)
    """

    u1 = r1/np.linalg.norm(r1)
    c23 = np.cross(r2, r3)
    c23 = c23/np.linalg.norm(c23)

    if abs(np.dot(u1, c23)) > tol:
        print('Cannot determine orbit: position vectors are not co-planar')
        return None

    N = np.linalg.norm(r1)*np.cross(r2, r3) + np.linalg.norm(r2)*np.cross(r3, r1) + np.linalg.norm(r3)*np.cross(r1, r2)
    D = np.cross(r1, r2) + np.cross(r2, r3) + np.cross(r3, r1)
    n = np.linalg.norm(N)
    d = np.linalg.norm(D)
    S = r1*(np.linalg.norm(r2) - np.linalg.norm(r3)) + r2*(np.linalg.norm(r3) - np.linalg.norm(r1)) \
        + r3*(np.linalg.norm(r1) - np.linalg.norm(r2))

    v2 = np.sqrt(mu/(n*d))*((np.cross(D, r2))/np.linalg.norm(r2) + S)

    return np.degrees(func.vector_to_elements(r2, v2, mu))


def lambert_bisection(r1, r2, Dt, mu):
    """ Function to solve Lambert's problem via bisection
        Based on bisection method in: http://control.asu.edu/Classes/MAE462/462Lecture10.pdf

        :param r1: initial position vector (m)
        :param r2: final position vector (m)
        :param Dt: time-of-flight (seconds)
    """

    # set up constants
    c = np.linalg.norm(r2 - r1)  # chord
    s = (c + np.linalg.norm(r1) + np.linalg.norm(r2)) / 2  # semi-perimeter
    amin = s / 2
    amax = 2 * s
    n = 0

    # check feasibility of solve
    tmin = np.sqrt(2) / 3 * np.sqrt(s ** 3 / mu) * (1 - ((s - c) / s) ** (3 / 2))

    if tmin > Dt:
        print('Cannot solve lambert\'s problem for Dt = ' + str(Dt) +
              ' s: TOF too small. Minimum time of: ' + str(tmin) + ' s required.')
        return None

    # iterate to solve for semimajor axis
    while True:
        a = (amax + amin) / 2
        alpha = np.arcsin(np.sqrt(s / (2 * a))) * 2
        beta = np.arcsin(np.sqrt((s - c) / (2 * a))) * 2

        g = np.sqrt(a ** 3 / mu) * (alpha - beta - (np.sin(alpha) - np.sin(beta)))

        if abs(g - Dt) <= 1e-6:
            break

        if n >= 1000:
            break

        if g > Dt:
            amin = a
        elif g < Dt:
            amax = a

        n += 1

    # calculate velocities
    A = np.sqrt(mu / (4 * a)) * (1 / np.tan(alpha / 2))
    B = np.sqrt(mu / (4 * a)) * (1 / np.tan(beta / 2))

    u1 = r1 / np.linalg.norm(r1)
    u2 = r2 / np.linalg.norm(r2)
    uc = (r2 - r1) / c

    v1 = (B + A) * uc + (B - A) * u1
    v2 = (B + A) * uc - (B - A) * u2

    return v1, v2


# # Earth-Mars transfer
# b1 = ephemerides.earth
# b2 = ephemerides.mars
#
# au = 1.496e11
# mu = 1.32712440018e20
#
# t1 = time.date_to_jd(30, 7, 2020, hour=12, mins=53)
# t2 = time.date_to_jd(30, 2, 2021, hour=13)
# ts = np.linspace(t1, t2, 2000)
# Dt = (t2 - t1)*24*60**2
#
# r1 = ephemerides.ephemeris(b1, t1)*au
# r2 = ephemerides.ephemeris(b2, t2)*au
#
# v = lambert_bisection(r1, r2, Dt, mu)
# v12, v22 = lambert(r1, r2, mu, Dt)
#
# # Plot the transfer in X-Y plane
# if v is not None:
#     v1, v2 = v
#
#     a, e, i, omega, w, M = func.vector_to_elements(r1, v1, mu)
#     t1 = func.orbit_from_elements(a, e, i, omega, w, M)
#     x, y, z = t1.plot()
#     plt.plot(np.array(x)/au, np.array(y)/au)
#
# # print(func.vector_to_elements(r1, v12, mu))
#
# a, e, i, omega, w, M = func.vector_to_elements(r1, v12, mu)
# print(func.orbit_from_elements(a, e, i, omega, w, M)
# x, y, z = t2.plot()
# plt.plot(np.array(x)/au, np.array(y)/au)
#
# # Plot planets and trajectory
# planets = [b1, b2]
#
# for T in ts:
#     for body in planets:
#         r = ephemerides.ephemeris(body, T)
#         body['x'].append(r[0])
#         body['y'].append(r[1])
#         body['z'].append(r[2])
#
# for p in planets:
#     plt.plot(p['x'], p['y'])
#
# plt.plot([0, r1[0]/au], [0, r1[1]/au], 'k')
# plt.plot([0, r2[0]/au], [0, r2[1]/au], 'k--')
# plt.axis('equal')
# plt.show()

