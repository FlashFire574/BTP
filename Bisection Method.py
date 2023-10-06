import numpy as np
from math import *

def func(x):
    num=x**2-9
    #num=np.exp(x-9)-20
    return num
epsilon=10**(-3)

a=float(input('enter guess root a:'))
b=float(input('enter guess root b:'))

while func(a)*func(b)>0:
    a=float(input('enter guess root a:'))
    b=float(input('enter guess root b:'))

i=0
while abs(a-b)>epsilon:
    c=(a+b)/2
    if func(a)*func(c)>0:
        a=c
    else:
        b=c
    i=i+1
root=(a+b)/2
print(root,i)
print(func(root))




