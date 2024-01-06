# Solving Poisson equation on a triangular domain with u(x) = 0 on the boundary

import numpy as np

def f(x,y):
    return -16*np.sin(y) + 8*(2*y-np.sqrt(3))*np.cos(y) - ((2*y-np.sqrt(3))**2)*np.sin(y) + 3*((2*x-1)**2)*np.sin(y)

def u(x,y):
    return ((2*y-np.sqrt(3))**2 - 3*((2*x-1)**2))*np.sin(y)

def getIndex(i,j,n):
    
    if j == 1:
        return i-1
    else:
        sum = 0
        for k in range (1,j):
            sum = sum + (n-2) - (k-1)
        return sum + (i-1)


def makeA(n):

    h = 1/n
    numPoints = int((n-2)*(n-1)/2)

    diag = np.ones(numPoints)*(-4/(h**2))
    A = np.diag(diag)

    for j in range(1,n-1):
        for i in range(1,n-j):

            if i+1 < n - j:
                A[getIndex(i,j,n)][getIndex(i+1,j,n)] = 4/(3*h**2) # right
            if j - 1 != 0:
                A[getIndex(i,j,n)][getIndex(i,j-1,n)] = 4/(3*h**2) # down left
            if i - 1 != 0:
                A[getIndex(i,j,n)][getIndex(i-1,j+1,n)] = 4/(3*h**2) # up left
        
    return A

def makeB(n):
    
    h = 1/n
    b = np.zeros(int((n-2)*(n-1)/2))

    for j in range(1,n-1):
        for i in range(1,n-j):
            x1 = h*(i +(1/2)*j)
            y = h*j*np.sqrt(3)/2
            b[getIndex(i,j,n)] = f(x1,y)

    return b

def error(x,n):
    
    error = 0
    for j in range(1,n-1):
        for i in range (1,n-j):
            x1 = h*(i+(1/2)*j)
            y = h*j*np.sqrt(3)/2        
            error = error + (x[getIndex(i,j,n)] - u(x1,y))**2

    error = np.sqrt((np.sqrt(3)/(4*n**2))*error)
    return error

nValues = [10,20,40,80,160]

for n in nValues:
    h = 1/n
    
    A = makeA(n)
    b = makeB(n)

    x = np.linalg.solve(A,b)
    print(n, error(x,n))
