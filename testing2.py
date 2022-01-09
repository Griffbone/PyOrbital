import numpy as np
import functions as func
import constants as cns
import matplotlib.pyplot as plt

# a, e, i, lan, w, ta,
rc, vc = func.elements_to_vector(6378e3, 0, 180, 0, 0, 135, cns.mu)

print(rc)
