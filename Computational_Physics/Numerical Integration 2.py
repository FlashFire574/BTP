from math import *

def f(x,y):
    z=x*(y**2)
    return z

sumy=0
a=2.
b=1.

n=400

dx=a/(n-1)
dy=b/(n-1)

y=0.+dy

for j in range(n-2):
    sumx=0
    x=0.+dx
    for i in range(n-2):
        sumx=sumx+f(x,y)
        x=x+dx
    sumx=(2*sumx+f(0,y)+f(a,y))*dx/2
    sumy=sumy+sumx
    y=y+dy
y=0
sumx=0
x=0.+dx

for i in range(n-2):
    sumx=sumx+f(x,y)
    x=x+dx
sumy0=(2*sumx+f(0,y)+f(a,y))*dx/2

y=b
sumx=0
x=0.+dx

for i in range(n-2):
    sumx=sumx+f(x,y)
    x=x+dx
sumyn=(2*sumx+f(0,y)+f(a,y))*dx/2

sumy=(2*sumy+sumy0+sumyn)*dy/2

print(sumy)










