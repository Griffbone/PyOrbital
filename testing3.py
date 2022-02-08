import numpy as np
import functions as func
import astrotime as at

# r = [0.36111 0.62563 -1.25093];		% DU
# rhat = r/norm(r);
# v = [0.37979 1.067295 0.263132];	% DU/TU

r = np.array([0.36111, 0.62563, -1.25093])
v = np.array([0.37979, 1.067295, 0.263132])

[a, e, i, lan, w, ta, arglat, truelon, lonper] = func.vector_to_elements(r, v, 1)

print(a)
print(e)
print(i)
print(lan)
print(w)
print(ta)
