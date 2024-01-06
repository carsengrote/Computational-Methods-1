
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

def makeMinv(n):

    numPoints = int((n-2)*(n-1)/2)

    diag = np.ones(numPoints)*(-4/(h**2))
    A = np.diag(diag)

    for j in range(1,n-1):
        for i in range(1,n-j):

            if i+1 < n - j:
                A[getIndex(i,j,n)][getIndex(i+1,j,n)] = 2/(3*h**2) # Right
            if j - 1 != 0:
                A[getIndex(i,j,n)][getIndex(i,j-1,n)] = 2/(3*h**2) # Down left
            if i - 1 != 0:
                A[getIndex(i,j,n)][getIndex(i-1,j+1,n)] = 2/(3*h**2) # up left
            if i - 1 != 0:
                A[getIndex(i,j,n)][getIndex(i-1,j,n)] = 2/(3*h**2) # left
            if i+1 < n - j: 
                A[getIndex(i,j,n)][getIndex(i,j+1,n)] = 2/(3*h**2) # up right
            if j -1 != 0: 
                A[getIndex(i,j,n)][getIndex(i+1,j-1,n)] = 2/(3*h**2) # down right

    blockSize = int(np.floor(np.sqrt(numPoints)))
    #blockSize = int(np.floor(numPoints/2))
    #fullInv = np.linalg.inv(A)
    
    lastBlockSize = numPoints % blockSize
    numBlocks = int(np.floor(numPoints/blockSize))
    Minv = np.zeros((numPoints,numPoints))
    #print(blockSize, lastBlockSize, numBlocks)
    #M = np.zeros((numPoints,numPoints))

    for i in range(0,numBlocks):
        Minv[i*blockSize:(i+1)*blockSize,i*blockSize:(i+1)*blockSize] = np.linalg.inv(A[i*blockSize:(i+1)*blockSize,i*blockSize:(i+1)*blockSize])
        #M[i*blockSize:(i+1)*blockSize,i*blockSize:(i+1)*blockSize] = A[i*blockSize:(i+1)*blockSize,i*blockSize:(i+1)*blockSize]
        

    if lastBlockSize > 0:
        Minv[numBlocks*blockSize:,numBlocks*blockSize:] = np.linalg.inv(A[numBlocks*blockSize:,numBlocks*blockSize:])
        #M[numBlocks*blockSize:,numBlocks*blockSize:] = A[numBlocks*blockSize:,numBlocks*blockSize:]
        
    return Minv

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
    h = 1/n
 
    mIv = makeMinv(n)

    u = np.zeros(int((n-2)*(n-1)/2))
    r = makeF(n)
    p = np.matmul(mIv,makeF(n))
    y = np.matmul(mIv,r)
    t = 0
    itCount = 0

    while(t==0):

        z = mulA(p,n)
        v = np.dot(y,r)/np.dot(p,z)
        u = u + v*p
        rNew = r - v*z

        #print(np.linalg.norm(rNew))
        if np.linalg.norm(rNew) < 10**(-10):
            t = 1
        
        yNew = np.matmul(mIv,rNew)
        mu = np.dot(yNew,rNew)/np.dot(y,r)
        p = yNew + mu*p
        r = rNew
        y = yNew

        itCount+=1

    print(n, itCount)

#print(np.linalg.norm(mulA(u,n) - makeF(n)))
#print(u)
