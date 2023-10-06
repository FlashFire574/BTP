from math import *
import numpy as np

def func(x):
    #num=x**2-9
    num=a*x**2+b*x+c
    return num
def grad_func(x):
    grad_num=2*a*x+b
    return grad_num
epsilon=10**(-3)
x=float(input('Enter x:'))
a=float(input('Enter a:'))
b=float(input('Enter b:'))
c=float(input('Enter c:'))

i=0
while abs(func(x))>epsilon:
    x=x-(func(x)/grad_func(x))
    i=i+1
root=x
print(root,i)
print(func(x))



