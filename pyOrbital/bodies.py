import numpy as np
import pyOrbital.constants as cn


class MajorBody:
    def __init__(self, re, rm, mass, rot):
        """Initialize a major body
            :param re: equitorial radius (m)
            :param rm: mean radius (m)
            :param mass: mass (kg)
            :param rot: sideral rotation period (days)
        """

        self.re = re
        self.rm = rm
        self.mass = mass
        self.rot = rot
        self.mu = cn.G*mass

    def vcircular(self, r):
        """Function to determine circular orbit velocity
            :param r: circular orbit radius (m)
        """

        return np.sqrt(self.mu/r)

    def vescape(self, r):
        """Function to determine escape velocity
            :param r: distance from body center (m)
        """

        return np.sqrt(2*self.mu/r)


# https://ssd.jpl.nasa.gov/?planet_phys_par
mercury = MajorBody(2440.53e3, 2439.4e3, 0.330144e24, 58.6462)
earth = MajorBody(6378.1366e3, 6371.0084e3, 5.97237e24, 0.99726968)

# print(earth.vescape(earth.rm))
