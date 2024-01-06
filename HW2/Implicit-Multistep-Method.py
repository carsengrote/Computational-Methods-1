
import numpy as np

def yExact(t):
   return np.exp(-t)*(np.cos(5*t) + (1/5)*np.sin(5*t))

def ypExact(t):
   return (-26/5)*np.exp(-t)*np.sin(5*t)

def yppExact(t):
   return (26/5)*np.exp(-t)*(np.sin(5*t)-5*np.cos(5*t))

hValues = [1/10,1/50,1/100,1/200,1/300,1/400,1/500,1/600,1/700,1/800,1/900,1/1000]

A = np.zeros((2,2))
A[0,0] = 0
A[0,1] = 1
A[1,0] = -26
A[1,1] = -2

for h in hValues:

    n = 3/h + 1
    n = int(n)
    if n % 10 == 0:
       n = n+1

    sol = np.zeros((2,n))
    sol[0,0] = 1
    sol[1,0] = 0
    #sol[0,1] = yExact(h)
    #sol[1,1] = ypExact(h)
    #sol[0,2] = yExact(2*h)
    #sol[1,2] = ypExact(2*h)
    
    sol[0,1] = sol[0,0]+h*ypExact(0)
    sol[1,1] = sol[1,0]+h*yppExact(0)
    sol[0,2] = sol[0,1]+h*ypExact(h)
    sol[1,2] = sol[1,1]+h*yppExact(h)

    inv = np.identity(2) - (h/3)*A
    inv = np.linalg.inv(inv)

    for k in range(2,n-1):
        
      step = np.matmul(inv,sol[:,k-1] + (h/3)*(np.matmul(A,sol[:,k-1]) + 4*np.matmul(A,sol[:,k])))
      sol[:,k+1] = step

    print(h,np.abs(yExact(3) - sol[0,-1]))
