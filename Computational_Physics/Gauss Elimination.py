from math import*
import numpy as np

#A=[[1,2,3],[4,5,6]] #string
#A=np.array(A)
#print(A*10)
#print(A[1][:])

A=np.array([[2.,3.,-1.],[4.,4.,-3.],[-2.,3.,-1.]])
A0=np.array([[2.,3.,-1.],[4.,4.,-3.],[-2.,3.,-1.]])
I=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
B=np.array([[5.],[3.],[1.]])
n=len(A)

for i in range(n):
    imx=np.argmax(A[i:n,i])
    mx=max(A[i:n,i])
    A[i+imx,:]/=mx
    I[i+imx,:]/=mx
    A[[i,i+imx],:]=A[[i+imx,i],:]
    I[[i,i+imx],:]=I[[i+imx,i],:]
    for j in range(i+1,n):
        I[j,:]-=A[[j],[i]]*I[i,:]
        A[j,:]-=A[[j],[i]]*A[i,:]
print(A0)
print(A)
print(I)
print(np.dot(I,B))
