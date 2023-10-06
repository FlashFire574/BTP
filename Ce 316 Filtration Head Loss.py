g=9.81
Pw=float(input('Density of water'))
dp=float(input('Diameter'))
s=float(input('Shape factor'))
u=float(input('Dynamic visocity'))
Vs=float(input('Filtering Velocity'))
e=float(input('Porosity'))
L=float(input('Depth of bed'))


Re= s*Vs*Pw*dp/u
print('Re = {}'.format(Re))

if Re<1:
    print("Laminar Flow Re = {}".format(Re))
    f=150*(1-e)/Re + 1.75
    print('f = {}'.format(f))
    hf=f*L*(1-e)*Vs*Vs/(e*e*e*g*dp)
    print('hf= {}'.format(hf))


