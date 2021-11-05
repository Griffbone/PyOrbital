import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import functions as func
import constants as cns
from matplotlib.patches import Circle
from matplotlib import animation
import matplotlib

import propagators

ra = 0.3744e9
rp = 200000 + cns.re
sat = func.Elements((ra + rp)/2, (ra - rp)/(ra + rp), 45, 0, 0, 0, cns.mu)
moon = func.Elements(0.3844e9, 0, 0, 0, 0, 140, cns.mu)

n = 1000
ts, xs, ys, _ = propagators.kepler_propagation(sat.a, sat.e, sat.i, 0, 0, sat.ta, moon.T, n)
_, xm, ym, _ = propagators.kepler_propagation(moon.a, moon.e, moon.i, 0, 0, moon.ta, moon.T, n)

########################################################################################################
fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(autoscale_on=False, xlim=(-moon.a*1.2, moon.a*1.2), ylim=(-moon.a*1.2, moon.a*1.2))
ax.set_aspect('equal')
ax.add_patch(Circle((0, 0), 6378e3, alpha=0.5))

########################################################################################################
dot, = plt.plot([xs[0]], [ys[0]], 'rx')
dot2, = plt.plot([xm[0]], [ym[0]], 'bx')

soi = [Circle((xm[0], ym[0]), 66.1e6, facecolor='w', edgecolor='k', linestyle='--')]
ax.add_patch(soi[0])

mxs = []
mys = []
mline, = plt.plot([], [], 'k--',  linewidth=0.25)

sxs = []
sys = []
sline, = plt.plot([], [], 'k--', linewidth=0.25)


def update(i):
    if i == 0:
        mxs.clear()
        mys.clear()
        sxs.clear()
        sys.clear()

    soi[0].remove()
    soi[0] = Circle((xm[i], ym[i]), 66.1e6, facecolor='w', edgecolor='k', linestyle='--')
    ax.add_patch(soi[0])

    dot.set_data([xs[i]], ys[i])
    dot2.set_data(xm[i], ym[i])

    mxs.append(xm[i])
    mys.append(ym[i])
    mline.set_data(mxs, mys)

    sxs.append(xs[i])
    sys.append(ys[i])
    sline.set_data(sxs, sys)

    return soi[0], dot, dot2, mline, sline


matplotlib.rcParams['animation.ffmpeg_path'] = r'C:\Users\Griffin\Downloads\ffmpeg-2021-11-03-git-08a501946f-essentials_build\ffmpeg-2021-11-03-git-08a501946f-essentials_build\bin\ffmpeg.exe'

ani = FuncAnimation(fig, update, n, interval=1, blit=True, repeat=False, save_count=50)
writerg = animation.FFMpegWriter(fps=100)
ani.save('test.mp4', writer=writerg)
plt.show()
