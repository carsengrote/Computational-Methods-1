import numpy as np
import sys


def f(x,y):
    return -16*np.sin(y) + 8*(2*y-np.sqrt(3))*np.cos(y) - ((2*y-np.sqrt(3))**2)*np.sin(y) + 3*((2*x-1)**2)*np.sin(y)


def getIndex(i,j,n):
    
    if j == 1:
        return i-1
    else:
        sum = 0
        for k in range (1,j):
            sum = sum + (n-2) - (k-1)
        return sum + (i-1)
    
def getIJ(k,n):
    i = 0
    j = 0
    rowLength = n - 2
    while(k/rowLength >= 1):
        j = j + 1
        k = k - rowLength
        rowLength = rowLength - 1
    i = k
    return i+1,j+1

def makeF(n):
    
    h = 1/n
    b = np.zeros(int((n-2)*(n-1)/2))

    for j in range(1,n-1):
        for i in range(1,n-j):
            x = h*(i+(1/2)*j)
            y = h*j*np.sqrt(3)/2
            b[getIndex(i,j,n)] = f(x,y)

    return b

def mulA(u,n):

    h = 1/n
    numPoints = int((n-2)*(n-1)/2)
    b = np.zeros(numPoints)
    c = 2/(3*h**2)
    
    for k in range(0,numPoints):

        b[k] = (-4/(h**2))*u[k]
        i,j = getIJ(k,n)

        if i+1 < n - j:
            b[k] +=  c*u[getIndex(i+1,j,n)] + c*u[getIndex(i,j+1,n)]
        if j - 1 != 0:
            b[k] += c*u[getIndex(i,j-1,n)] + c*u[getIndex(i+1,j-1,n)]
        if i - 1 != 0:
            b[k] += c*u[getIndex(i-1,j,n)] + c*u[getIndex(i-1,j+1,n)]

    return b



nValues = (sys.argv[1]).split()

for n in nValues:
    n = int(n)
    u = np.zeros(int((n-2)*(n-1)/2))
    r = makeF(n) - mulA(u,n)
 
    p = r
    t = 0
    itCount = 0
    while(t==0):
 
        w = mulA(p,n)

        alpha = np.dot(r,r)/np.dot(p,w)

        u = u + alpha*p
        rNew = r - alpha*w

        if np.linalg.norm(rNew) < 10**(-10):
            t = 1
        beta = np.dot(rNew,rNew)/np.dot(r,r)
        p = rNew + beta*p
        r = rNew
        itCount+=1

    print(n, itCount)

#print(np.linalg.norm(mulA(u,n) - makeF(n)))
#print(u)