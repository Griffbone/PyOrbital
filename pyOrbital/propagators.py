import numpy as np
from scipy.optimize import fsolve


def kepler_E(e, M):
    """Solve Kepler's equation for Eccentric anomaly
        :param e: eccentricity
        :param M: mean anomaly
        :return E: eccentric anomaly
    """

    if (e < 0) or (e > 1):
        print('ERROR: eccentricty must be below one')
        return None
    elif (M < 0) or (M > np.pi):
        print('ERROR: mean anomaly must be below one')
        return None

    return fsolve(lambda E: E - e*np.sin(E) - M, np.array([0]))


# def kepler_propagate(DT):


# e = 0.5
# M = np.radians(40)
# print(np.degrees(kepler_E(e, M)))
# import constants as cn
# import numpy as np
# from scipy.optimize import fsolve
#
# global CK2, XJ3, AE, XKE
# CK2 = 5.413080E-4
# XJ3 = -.253881E-5
# AE = 1.0
# XKE = .743669161E-1
#
#
# def FMOD2P(X):
#     FMOD=X
#     I=FMOD/(np.pi*2)
#     FMOD=FMOD-I*np.pi*2
#
#     if FMOD < 0:
#         FMOD = FMOD*np.pi*2
#
#     return FMOD
#
#
# def sgp(XINCL, XN0, E0, XM0, OMEGA0, XNODEO, XNDT20, XNDD60, TSINCE):
#     """Function to propogate TLEs with SGP method
#         :param XINCL: mean inclination at epoch
#         :param XNO: mean motion at epoch
#         :param E0: mean eccentricity at epoch
#         :param XMO: mean anomaly at epoch
#         :param OMEGA0: mean argPe at epoch
#         :param XNODEO: mean RAAN at epoch
#         :param XND: first derivative of mean motion
#         :param XNDD: second derivative of mean motion
#         :param TSINCE: time since epoch
#     """
#
#     XNDT20 = XNDT20/2
#     XNDD60 = XNDD60/6
#
#     # constants
#     C1 = CK2*1.5
#     C2 = CK2/4
#     C3 = CK2/2
#     C4 = XJ3*AE**3/(4*CK2)
#
#     COSI0 = np.cos(XINCL)
#     SINI0 = np.sin(XINCL)
#     A1 = (XKE/XN0)**(2/3)
#     D1 = C1/(A1**2)*(3*COSI0**2 - 1)/(1 - E0**2)**(3/2)
#     A0 = A1*(1 - (1/3)*D1 - D1**2 - (134/81)*D1**3)
#     P0 = A0*(1 - E0**2)
#     Q0 = A0*(1 - E0)
#     XLO = XM0 + OMEGA0 + XNODEO
#
#     D10 = C3*SINI0**2
#     D20 = C2*(7*COSI0**2 - 1)
#     D30 = C1*COSI0
#     D40 = D30*SINI0
#
#     P02N0 = XN0/(P0**2)
#     OMGDT = C1*P02N0*(5*COSI0**2 - 1)
#     XNODOT = -2*D30*P02N0
#     C5 = (1/2)*C4*SINI0*(3 + 5*COSI0)/(1 + COSI0)
#     C6 = C4*SINI0
#
#     # update for gravity and drag
#     A = XN0 + (2*XNDT20 + 3*XNDD60*TSINCE)*TSINCE
#     A = A0*(XN0/A)**(2/3)
#
#     if A > Q0:
#         E = 1 - Q0/A
#     else:
#         E = 1e-6
#
#     P = A*(1 - E**2)
#     XNODES = XNODEO + XNODOT*TSINCE
#     OMGAS = OMEGA0 + OMGDT*TSINCE
#     XLS = FMOD2P(XLO + (XN0 + OMGDT + XNODOT + (XNDT20 + XNDD60*TSINCE)*TSINCE)*TSINCE)
#     AXNSL = E*np.cos(OMGAS)
#     AYNSL = E*np.sin(OMGAS) - C6/P
#     XL = FMOD2P(XLS-C5/P*AXNSL)
#
#     # solve kepler's equation
#     U = FMOD2P(XL - XNODES)
#     ITEM3 = 0
#     E01 = U
#     TEM5 = 1
#
#     while True:
#         SINE01 = np.sin(E01)
#         COSE01 = np.cos(E01)
#
#         if abs(TEM5) < 1e-6:
#             break
#         elif ITEM3 > 10:
#             break
#
#         ITEM3 += 1
#         TEM5 = 1 - COSE01*AXNSL - SINE01*AYNSL
#         TEM5 = (U - AYNSL*COSE01 + AXNSL*SINE01 - E01)/TEM5
#         TEM2 = abs(TEM5)
#
#         if TEM2 > 1:
#             TEM5 = TEM2/TEM5
#
#         E01 = E01 + TEM5
#
#     # short period prelim quantities
#     ECOSE = AXNSL*COSE01 + AYNSL*SINE01
#     ESINE = AXNSL*SINE01 - AYNSL*COSE01
#     EL2 = AXNSL**2 + AYNSL**2
#     PL = A*(1 - EL2)
#     PL2 = PL**2
#     R = A*(1 - ECOSE)
#     RDOT = XKE*np.sqrt(A)/R*ESINE
#     RVDOT = XKE*np.sqrt(PL)/R
#     TEMP = ESINE/(1 + np.sqrt(1 - EL2))
#     SINU = A/R*(SINE01 - AYNSL - AXNSL*TEMP)
#     COSU = A/R*(COSE01 - AXNSL + AYNSL*TEMP)
#     SU = np.arctan2(SINU, COSU)                 # ACTAN in the FORTRAN implementation
#
#     # short periodics
#     SIN2U = COSU*2*SINU
#     COS2U = 1 - 2*SINU**2
#     RK = R + D10/PL*COS2U
#     UK = SU - D20/PL2*SIN2U
#     XNODEK = SU - D20/PL2*SIN2U
#     XINCK = XINCL + D40/PL2*COS2U
#
#     # orientation vectors
#     SINUK = np.sin(UK)
#     COSUK = np.cos(UK)
#     SINNOK = np.sin(XNODEK)
#     COSNOK = np.cos(XNODEK)
#     SINIK = np.sin(XINCK)
#     COSIK = np.cos(XINCK)
#     XMX = -SINNOK*COSIK
#     XMY = COSNOK*COSIK
#
#     UX = XMX*SINUK + COSNOK*COSUK
#     UY = XMY*SINUK + SINNOK*COSUK
#     UZ = SINIK*SINUK
#
#     VX = XMX*COSUK - COSNOK*SINUK
#     VY = XMY*COSUK - SINNOK*SINUK
#     VZ = SINIK*COSUK
#
#     # position and velocity
#     X = RK*UX
#     Y = RK*UY
#     Z = RK*UZ
#
#     XDOT = RDOT*UX
#     YDOT = RDOT*UY
#     ZDOT = RDOT*UZ
#
#     XDOT = RVDOT*VX + XDOT
#     YDOT = RVDOT*VY + YDOT
#     ZDOT = RVDOT*VZ + ZDOT
#
#     return np.array([X, Y, Z]), np.array([XDOT, YDOT, ZDOT])
#
#
# XMNPDA = 1440
# TEMP = 2*np.pi/XMNPDA**2
# XKMPER = 6378.135
# E0 = 0.0086731                  # eccentricity
#
#
# # Incl=720435 RghtAscOfAscNode=15.689 Ecc=086731 ArguPerigee=206880 MeanAnomaly=005140 MeanMotion=00824518
# XINCL = np.radians(72.0435)         # inclination
# XNODEO = np.radians(15.689)         # RAAN
# OMEGA0 = np.radians(52.6988)        # ArgPe
# XM0 = np.radians(110.5714)          # Mean anomaly
# XNDT20 = .00073094*TEMP             # first derivative of mean motion
# XNDD60 = 13844e-3*TEMP/XMNPDA       # second derivative of mean motion
# XN0 = 16.05824518*TEMP*XMNPDA       # mean motion
#
# TSINCE = 0
# Rvec, Vvec = sgp(XINCL, XN0, E0, XM0, OMEGA0, XNODEO, XNDT20, XNDD60, TSINCE)
#
# Rvec = Rvec * XKMPER/AE
# Vvec = Vvec * XKMPER/AE*XMNPDA/86400
#
# print(Rvec)
# print(Vvec)
#
# """Function to propogate TLEs with SGP method
#     :param XINCL: mean inclination at epoch
#     :param XNO: mean motion at epoch
#     :param E0: mean eccentricity at epoch
#     :param XMO: mean anomaly at epoch
#     :param OMEGA0: mean argPe at epoch
#     :param XNODEO: mean RAAN at epoch
#     :param XND: first derivative of mean motion
#     :param XNDD: second derivative of mean motion
#     :param TSINCE: time since epoch
# """
#
