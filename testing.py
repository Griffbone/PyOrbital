import numpy as np
import matplotlib.pyplot as plt
import constants as cns
import functions as func
import dev_functions as dev


def circle(r, cx, cy, style):
    ts = np.linspace(0, 2*np.pi, 100)
    xs = r*np.cos(ts) + cx
    ys = r*np.sin(ts) + cy

    plt.plot(xs, ys, style)


h = 80000*1000**2
e = 1.4
a = h**2/(cns.mu*(e**2 - 1))

a1 = np.arccos(-1/e)
print(a1)

_, _, x, y = func.perifocal_coords(a, e)
plt.plot(x, y)

x, y, _ = func.elements_to_eci_pos(a, e, 0, 0, 0, 90)
vx, vy, _ = func.elements_to_eci_vel(a, e, 0, 0, 0, 90)

plt.scatter(x, y)
plt.plot([x, x + vx*1000], [y, y + vy*1000])

plt.axis('equal')
plt.show()
