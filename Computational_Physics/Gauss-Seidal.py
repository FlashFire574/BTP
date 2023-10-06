from math import *
import numpy as np

A=np.array([[8., 3., -3.], [-2., -8., 5.], [3., 5., 10.]])
(x1,x2,x3)=0,0,0
X0=np.array([x1,x2,x3])
B=np.array([14.,5.,-8.])
epsilon=0.01

for k in range(50):
    x1= (14-3*x2+3*x3)/8
    x2=(5+2*x1-5*x3)/(-8)
    x3 = (-8-3*x1-5*x2)/(-5)
    X = np.array([x1, x2, x3])
    print(k,X)
    if np.sqrt(np.dot(X-X0,X-X0))< epsilon:
        break
    else:
        X0=X








