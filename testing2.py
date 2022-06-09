import astrotime as ast
import functions as func
import numpy as np

year = 2022
month = 2
day = 11    # day after the 10th
hour = 0    # 17:00 MST = 24:00 UTC
min = 35

[jdn, jdf] = ast.date_to_jd(year, month, day, hour, min, 0)

print(jdn)
print(jdf)
