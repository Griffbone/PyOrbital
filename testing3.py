import numpy as np
import astrotime as ast

jdn, jdf = ast.date_to_jd(1992, 8, 20, 12, 14, 0)

gmst = np.rad2deg(ast.jd_to_gmst(jdn + jdf))
lst = np.rad2deg((ast.jd_to_lst(jdn + jdf, -104)))

print(ast.deg_to_hms(gmst))
print(ast.deg_to_hms(lst))
