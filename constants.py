# (1) https://ssd.jpl.nasa.gov/astro_par.html
# (2) https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf
# (3) https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
# (4) https://confluence.qps.nl/qinsy/latest/en/world-geodetic-system-1984-wgs84-182618391.html

G = 6.67430e-11         # universal gravitational constant - m**3/kg/s**2       (1)
au = 149597870700       # one AU - m                                            (1)
c = 299792458           # speed of light - m/s                                  (1)
d = 86400               # day - s                                               (1)
j2000 = 2451545.0       # J2000 epoch julian date - days                        ( )

g0 = 9.80665            # Earth standard gravity - m/s**2   ( )
mu = 3.98600435507e14   # Earth mu - m**3/s**2              (2)
re = 6378.137e3         # Earth radius - m                  (3)
smae = 6378137.0        # Earth ellipse semimajor axis - m  (4)
flat = 298.257223563    # Earth flattening (1/f)            (4)
rot = 7292115e-11       # Earth angular velocity - rad/s    (4)
j2 = 1.08262668e-3      # Earth J2 constant                 ( )

mu_sun = 1.32712440041279419e20     # Sun mu - m**3/s**2    (2)
