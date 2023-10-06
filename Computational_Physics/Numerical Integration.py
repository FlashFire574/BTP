from math import *

def func(x):
    num=sqrt(1-(sin(x)*sin(x)))
    return num

sum=0
a=float(input('Enter a'))
b=float(input('Enter b'))
x=a
n=3
dx=(b-a)/(n-1)
for i in range(n):
     if x==a or x==b:
       sum= sum +func(x)*dx/3
     elif i%2==0:
       sum= sum +4*func(x)*dx/3
     elif i%2!=0:
         sum= sum + 2*func(x)*dx/3
     x=x+dx
print(sum)
