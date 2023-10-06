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
c=10**5
i=0
while abs(func(c))>epsilon:
    c=b-func(b)*(b-a)/(func(b)-func(a))
    if func(a)*func(c)>0:
        a=c
    else:
        b=c
    i=i+1
root=c
print(root,i)
print(func(root))
