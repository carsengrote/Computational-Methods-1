import numpy as np

# Function written by Chris Rycroft
# Computes the (n+1) Chebyshev-Lobatto points at the associate derivative matrix
def cheb(n):

    # Chebyshev-Lobatto points
    x=np.cos(np.linspace(0,np.pi,n+1))

    # Weighting array
    c=np.ones((n+1))
    c[0]=2
    c[n]=2

    # Set up diagonal terms
    D=np.empty((n+1,n+1))
    D[0,0]=(2*n*n+1)/6
    for i in range(1,n):
        D[i,i]=-0.5*x[i]/(1-x[i]*x[i])
    D[n,n]=-D[0,0]

    # Set up off-diagonal terms
    for i in range(n+1):
        for j in range(n+1):
            if i!=j:
                D[i,j]=c[i]*(-1)**(i+j)/(c[j]*(x[i]-x[j]))
    return (x,D)

def u(x):
    return np.exp(x)*(x**2 - 1)

def f(x):
    return np.exp(5*x)*((x**2 -1)**5) + np.exp(x)*(x**2 + 4*x + 1)

def Pn(x, colPoints, uk):
    
    ans = 0
    index = 0
    for n in colPoints:
        
        numSum = 1
        denSum = 1
        for j in colPoints:
            if (j != n):
                numSum = numSum * (x - j)
                denSum = denSum * (n - j)

        ans = ans + uk[index]*numSum/denSum
    
        index = index + 1
    return ans

def g(n,uk):
 
    points, D1 = cheb(n)
    D2 = np.matmul(D1,D1)
    g = np.matmul(D2,uk)
    
    g[0] = 0
    g[n] = 0

    for i in range(1,n):
        g[i] = g[i] + uk[i]**5 - f(points[i])

    return g

def gPrime(n,uk):
    
    points, D1 = cheb(n)
    D2 = np.matmul(D1,D1)
     
    for i in range(0,n+1):
        D2[0][i] = 0
        D2[n][i] = 0
        D2[i][i] = D2[i][i] + 5*(uk[i]**4)
    D2[0][0] = 1
    D2[n][n] = 1
    return D2


def l2Error(xValues, colPoints, uk, dx):

    error = 0 
    for x in xValues:
        error = error + ((Pn(x,colPoints,uk) - u(x))**2)*dx

    error = error - (dx/2)*((Pn(-1, colPoints, uk) - u(-1))**2 + (Pn(1, colPoints, uk)- u(1))**2)
    error = np.sqrt(error)

    return error 

def error(uk, points,N):
    s = 0
    for i in range(0,N+1):
        s = s + (uk[i] - u(points[i]))**2
    
    return np.sqrt(s*(2/len(points))) # 2/len(points) is average h?

nValues = [4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64]
xValues = np.linspace(-1,1,200)
dx = 2/200

for N in nValues:
    
    points,D= cheb(N)
    uk = np.ones(N+1)
    uk[0] = 0
    uk[N] = 0
    # Initial guess
    for i in range(1,N):
        uk[i] = ((points[i]+1)**2)*(points[i]-1)

    deltaU = np.linalg.solve(gPrime(N,uk),g(N,uk))
    uk = uk - deltaU
    itCount = 1

    while error(uk,points,N) > 10**(-12) and itCount < 1000:
        deltaU = np.linalg.solve(gPrime(N,uk),g(N,uk))
        uk = uk - deltaU 
        itCount += 1
        #print(error(uk,points,N))
        #print(np.linalg.norm(g(N,uk)))

    print("# iterations: ", itCount)
    # for x in xValues:
    #     print(x,u(x),Pn(x, points,uk))
    #for i in range(0,N+1):
    #    print(points[i],u(points[i]),uk[i])
    #print(uk) 
    print(N, l2Error(xValues, points, uk, dx))

