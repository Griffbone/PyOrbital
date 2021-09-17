import nump as np


a1 = (ke/no)**(2/3)

del1 = (3/2)*(k2/a1**2)*((3*np.cos(io)**2 - 1)/(1 - eo**2)**(3/2))
ao = a1*(1 - (1/3)*del1 - del1**2 - (134/81)*del1**3)
po = ao*(1 - eo**2)
qo = ao*(1 - eo)
Lo = Mo + wo + omegao
domegadt = -(3/2)*J2*(ae**2/po**2)*no*np.cos(io)
dwdt = (3/4)*J2*(ae**2/po**2)*no*(5*np.cos(io)**2 - 1)
a = ao*(no/(no + 2*(nodot/2)*(t - to) + 3*(noddot/6)*(t - to)**2))**(2/3)

if a > qo:
    e = 1 - qo/a
else:
    e = 10e-6

p = a*(1 - e**2)
omegaso = omegao + domegadt*(t - to)
wso = wo + dwdt*(t - to)
Ls = Lo + (no + dwdt + domegadt)*(t - to) + (nodot/2)*(t - to)**2 + (noddot/6)*(t - to)**3

aynsl = e*np.sin(wso) - (1/2)*(J3/J2)*(ae/p)*np.sin(io)
axnsl = e*np.cos(wo)
L = Ls - (1/4)*(J3/J2)*(ae/p)*axnsl*np.sin(io)*((3 + 5*np.cos(io))/(1 + np.cos(io)))

iterations = 0
while iterations <= 25:
        


delo = (3/2)*(k2/ao**2)*((3*np.cos(io)**2 - 1)/(1 - eo**2)**(3/2))
nopp = no/(1 + delo)
aopp = ao/(1 - delo)

