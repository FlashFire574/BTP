import numpy as np
import matplotlib.pyplot as plt

# Question 1
Fck=25.
Fy=500.
Pt=[]
l=[]

for pt in range(20,81,1):
    pt=pt/100
    Pt.append(pt)

for i in range(600,1850,50):
    l.append(i)

rows,columns=(25,61)
D=[[0 for i in range(columns)] for j in range(rows)]

i=0
for K1 in range(600,1850,50): # K1=Mu/b
    j=0
    for pt in range(20,81,1):
        pt=pt/100
        K2=0.87*Fy*pt/100*(1-(Fy*pt)/(Fck*100))*(10**3)
        d=np.sqrt(K1/K2)
        D[i][j]=d
        j=j+1
    i=i+1

for i in range(25):
    plt.plot(Pt,D[i],label=l[i])
    plt.title('M25, Fe500')
    plt.xlabel('pt (%)')
    plt.ylabel('d (meter)')
plt.legend()
plt.show()

# Question 2
Fck=30.
Fy=415.
Pt=[]
l=[]

for pt in range(20,81,1):
    pt=pt/100
    Pt.append(pt)

for i in range(600,1850,50):
    l.append(i)

rows,columns=(25,61)
D=[[0 for i in range(columns)] for j in range(rows)]

i=0
for K1 in range(600,1850,50):
    j=0
    for pt in range(20,81,1):
        pt=pt/100
        K2=0.87*Fy*pt/100*(1-(Fy*pt)/(Fck*100))*(10**3)
        d=np.sqrt(K1/K2)
        D[i][j]=d
        j=j+1
    i=i+1

for i in range(25):
    plt.plot(Pt,D[i],label=l[i])
    plt.title('M30, Fe415')
    plt.xlabel('pt (%)')
    plt.ylabel('d (meter)')
plt.legend()
plt.show()



