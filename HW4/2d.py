import numpy as np
# Interval is [0,2pi) with periodic BCs
N = 40  # Grid points
h = 2*np.pi/N  # Grid spacing
v = 1/100  # Wave speed?

# Initial conditions u = exp(2sin(x))
u = np.zeros(N)
for i in range(0,N):
    u[i] = np.exp(2*np.sin(i*h))

A = np.zeros((N,N)) # Modified Beam Warming, dropping the second v^2 term
                    # To update u, just need to multiply it by a matrix

for i in range(0,N):
    A[i,i] = 1 - 3*(v/2)

    if (i == 0): # Wrap around
        A[i,-1] = 4*(v/2)
        A[i,-2] = -(v/2)
    elif (i == 1): # Wrap around
        A[i,i-1] = 4*(v/2)
        A[i,-1] = -(v/2)
    else: # Normal
        A[i,i-1] = 4*(v/2)
        A[i,i-2] = -(v/2)

def calcR(u,length):
    sum = 0
    for i in range (0,length):
        sum  = sum + u[i]**2
    sum = sum*(1/length)
    return np.sqrt(sum)

n = 10000 # Number of time steps

R = np.zeros(n+1) # Root mean squared value of the solution
R[0] = calcR(u,N)

for k in range(0,n):
    u = np.matmul(A,u)
    R[k+1] = calcR(u,N)

#for i in range(0,N):
    #print(i*h, u[i])

for i in range(0,n+1):
	print(i,R[i])
