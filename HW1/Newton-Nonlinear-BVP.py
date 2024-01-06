# Newtons method to solve nonlinear BVP u''(x) - 80cos(u(x)) = 0
# with u(0) = 0, u(1) = 10

import numpy as np

n = 100
h = 1/n

def f(xk):

    F = np.zeros(n-1)
    F[0] = (1/(h**2))*(xk[1]-2*xk[0]) - 80*np.cos(xk[0])
    
    for i in range(1,n-2):
       F[i] = (1/(h**2))*(xk[i+1]-2*xk[i]+xk[i-1]) - 80*np.cos(xk[i])

    F[n-2] = (1/(h**2))*(10-2*xk[n-2]+xk[n-3]) - 80*np.cos(xk[n-2])

    return F

def fprime(xk):
    
    x = np.zeros((n-1,n-1))

    x[0][0] = (-2/(h**2) + 80*np.sin(xk[0]))
    x[0][1] = 1/(h**2)
    x[n-2][n-3] = 1/(h**2)
    x[n-2][n-2] = (-2/(h**2) + 80*np.sin(xk[n-2]))

    for i in range(1,n-2):  
        x[i][i-1] = 1/(h**2)
        x[i][i] = (-2/(h**2) + 80*np.sin(xk[i]))
        x[i][i+1] = 1/(h**2)

    return x

xk = np.zeros(n-1)
deltax = np.linalg.solve(fprime(xk),f(xk))
xk = xk - deltax
itCount = 1

while(np.linalg.norm(deltax)/np.linalg.norm(xk) > (10**-10)):
    
    deltax = np.linalg.solve(fprime(xk),f(xk))
    xk = xk - deltax
    itCount += 1
    print(np.linalg.norm(f(xk)))

print(0,0)
for i in range(0,n-1):
    print(i*h + h,xk[i])

print(1,10)

print(xk[50],xk[51],xk[52])
print(itCount)
