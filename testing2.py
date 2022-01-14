import numpy as np
import functions as func
import constants as cns


print(np.isnan(np.nan))
# mu = cns.mu/1000**3
# re = 6378
# c45 = np.cos(np.radians(45))
#
# r = re*np.array([c45, c45, 0])
# v = 1.2*np.sqrt(mu/re)*np.array([-c45, c45, 0])
# a, e, i, lan, w, ta, arglat, truelon, lonper = func.vector_to_elements(r, v, mu)
#
# print('a    : {}'.format(a))
# print('e    : {}'.format(e))
# print('i    : {}'.format(np.radians(i)))
# print('LAN  : {}'.format(np.radians(lan)))
# print('w    : {}'.format(np.radians(w)))
# print('ta   : {}\n'.format(np.radians(ta)))
#
# print('arglat   : {}'.format(np.radians(arglat)))
# print('truelon  : {}'.format(np.radians(truelon)))
# print('lonper   : {}'.format(np.radians(lonper)))
