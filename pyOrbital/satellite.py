import numpy as np


class Satellite:
    def __init__(self, parent, a=None, e=0, i=0, lan=0, argPe=0, ta=0):
        """Function to initialize a satellite object
            :param parent: parent object
            :type parent: pyOrbital.bodies.MajorBody
            :param elements: dictionary of orbital elements
        """

        self.parent = parent

        if a is None:
            a = self.parent.rm

        self.a = a              # semimajor axis (m)
        self.e = e              # eccentricity
        self.i = i              # inclination (deg)
        self.lan = lan          # longitude of the ascending node (deg)
        self.argPe = argPe      # argument of periapsis (deg)
        self.ta = ta            # true anomaly (deg)

        self.T = np.sqrt(((4*np.pi**2)/parent.mu)*a**3)     # period (s)
        self.n = 360/self.T                                 # mean motion (deg/s)

    def propagate(self, method='kepler', orbs = 1):
        """Function to propagate the satellite by various methods
        """

        pass

    def plot_orbit(self, plot_type='polar'):
        """Function to plot a single non-perturbed orbit of the satellite
        """

        pass

