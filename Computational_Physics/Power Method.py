from math import *
import numpy as np

def normalize(Z):
    Z_max=max(abs(Z))
    Zn=Z/max(Z)
    return Z_max,Zn

X = np.array([1, 1,1])
A = np.array([[3,1,4],
              [0,2, 6],[0,0,5]])

for i in range(10):
    X=np.dot(A,X)
    l,X=normalize(X)
    print(l,X)


