# BVP u'' + u = f(x) on domain [0,pi] with mixed boundary conditions
# u'(0) - u(0) = 0 and u'(pi) + u(pi) = 0
# Second order accurate

import numpy as np

def u(x):
    return (1/2)*(-1*np.exp(x) - np.exp(np.pi)*np.cos(x) - np.exp(np.pi)*np.sin(x))

def e(x,n):
    
    h = np.pi/n
    sum = 0

    for j in range(0,n+1):
        if (j == 0 or j == n): 
            q = 1/2
        else:
            q = 1
            
        sum = sum + q*((x[j] - u(h*j))**2)

    sum = sum*h
    return np.sqrt(sum)

meshWidth = [20,40,80,160]

for n in meshWidth:
    
    h = np.pi/n
    b = np.zeros(n+1)

    A = np.zeros((n+1,n+1))
    for i in range (1,n):
        A[i][i-1] = 1
        A[i][i] = (h**2 - 2)
        A[i][i+1] = 1

        b[i] = -1*np.exp(i*h) 

    A = (1/h**2)*A

    A[0][0] = (-1*3/2 - h)*(1/h)
    A[0][1] = (2)*(1/h)
    A[0][2] = (-1/2)*(1/h)
    
    A[n][n-2] = (1/2)*(1/h)
    A[n][n-1] = (-2)*(1/h)
    A[n][n] = (3/2 + h)*(1/h)

    x = np.linalg.solve(A,b)   
    
    print(n,e(x,n))
