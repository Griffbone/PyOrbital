import maneuvers as man
import pyOrbital.bodies
import pyOrbital.constants as cns
import numpy as np
import matplotlib.pyplot as plt
import plotting as oplot
import matplotlib.patches as patches


class Orbit:
    def __init__(self, a, e, i, lan, w, ta=0, parent=pyOrbital.bodies.earth):
        self.a = a
        self.e = e
        self.i = i
        self.lan = lan
        self.w = w
        self.ta = ta

    def plot_xy(self, ax, tas=np.linspace(0, 2*np.pi, 100)):
        x, y, _ = oplot.plot_orbit(self.a, self.e, self.i, self.lan, self.w, tas)
        ax.plot(x, y)

    def plot_yz(self, ax, tas=np.linspace(0, 2*np.pi, 100)):
        _, y, z = oplot.plot_orbit(self.a, self.e, self.i, self.lan, self.w, tas)
        ax.plot(y, z)

    def plot_xz(self, ax, tas=np.linspace(0, 2*np.pi, 100)):
        x, _, z = oplot.plot_orbit(self.a, self.e, self.i, self.lan, self.w, tas)
        ax.plot(x, z)


T = 60*60*24

r1 = 6578e3
r2 = (cns.mu*T**2/((2*np.pi)**2))**(1/3)
a = (r1 + r2)/2
e = (r2 - r1)/(r1 + r2)

_, v1, _, v2 = man.vis_viva(r1, r2, cns.mu)

dv1 = v1 - np.sqrt(cns.mu/r1)
dv2 = man.plane_change(v2, np.radians(28), np.sqrt(cns.mu/r2))

initial = Orbit(r1, 0, np.radians(28), 0, 0)
transfer = Orbit(a, e, np.radians(28), 0, 0)
final = Orbit(r2, 0, 0, 0, 0)

plt.style.use('seaborn')
fig, ax = plt.subplots()

initial.plot_xy(ax)
transfer.plot_xy(ax)
final.plot_xy(ax)

ax.set_aspect('equal')
plt.show()




# plt.show()

# earth = patches.Circle((0, 0), cns.re, alpha=0.5)
# ax.add_artist(earth)

# x, y, z = oplot.plot_orbit(r1, 0, 28, 0, 0)
# plt.plot(x, y, z, '--')
# max1 = max([max(x), max(y), max(z)])
#
# x, y, z = oplot.plot_orbit(a, e, 28, 0, 0, np.linspace(0, np.pi, 100))
# plt.plot(x, y, z, '--')
# max2 = max([max(x), max(y), max(z)])
#
# x, y, z = oplot.plot_orbit(r2, 0, 0, 0, 0)
# plt.plot(x, y, z, '--')
# max3 = max([max(x), max(y), max(z)])
#
# gmax = max([max1, max2, max3])
#
# ax.scatter(r1, 0, 0, zorder=3)
# ax.scatter(-r2, 0, 0, zorder=3)
#
# ax.set_xlim([-gmax, gmax])
# ax.set_ylim([-gmax, gmax])
# ax.set_zlim([-gmax, gmax])
#
#
# u = np.linspace(0, 2*np.pi, 20)
# v = np.linspace(0, np.pi, 20)
#
# x = cns.re * np.outer(np.cos(u), np.sin(v))
# y = cns.re * np.outer(np.sin(u), np.sin(v))
# z = cns.re * np.outer(np.ones(np.size(u)), np.cos(v))
#
# ax.plot_surface(x, y, z, alpha=1, edgecolor='k', facecolor='b')
#
#
# ax.set_xlabel('ECI X [m]')
# ax.set_ylabel('ECI Y [m]')
# ax.set_zlabel('ECI Z [m]')
#
# ax.legend(['Initial Orbit', 'Transfer Orbit', 'Final Orbit',
#            'First Burn, DV = {:.0f} m/s'.format(dv1), 'Second Burn, DV = {:.0f} m/s'.format(dv2)])
# plt.title('GTO-1800 Transfer Orbit')
#
# # plt.axis('equal')
# plt.show()
