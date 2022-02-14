import astrotime as ast

jdn, ut = ast.date_to_jd(2016, 12, 0, 0, 0, 0)

print(ast.gmst(jdn))



