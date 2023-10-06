import numpy as np

g=9.81
Pp=float(input('Specific gravity'))
Pw=float(input('Density of water'))
s=float(input('Shape factor'))
u=float(input('Dynamic visocity'))
d=float(input('Diameter'))

Vt=(Pp-1)*g*Pw*d*d/(18*u)
print('Vt ={} m/s'.format(Vt))

Re= s*Vt*Pw*d/u
print('Re = {}'.format(Re))

if Re<1:
 print("Laminar Flow Re = {}".format(Re))

elif 1<Re<10**4:
 print('Transitional flow')
 i=0
 while i<3:
   Cd=24/Re +3/np.sqrt(Re)+0.34
   print('Cd = {}'.format(Cd))
   Vt= np.sqrt(4*(Pp-1)*g*d/(3*Cd))
   print('Vt ={} m/s'.format(Vt))
   Re= s*Vt*Pw*d/u
   print('Re = {}'.format(Re))
   i=i+1

elif Re>10**4:
 print('Turbulent Flow')
 Cd=0.4
 Vt= np.sqrt(4*(Pp-1)*g*d*d/(3*Cd))
 print('Vt ={} m/s'.format(Vt))





