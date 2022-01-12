import numpy as np
import functions as func
import constants as cns
import matplotlib.pyplot as plt

r = 6378*np.array([1, 0, 0])
v = 9*np.array([0, 1, 0])

a, e, i, lan, w, ta, truelon, arglat, lonper = func.vector_to_elements(r, v, cns.mu/1000**3)

print([a, e, i, lan, w, ta])
print([truelon, arglat, lonper])
