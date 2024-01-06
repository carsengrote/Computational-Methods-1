# Solving 2D diffusion with Alternating Direction Implicit method
# With Dirichilet BCs and IC u(x, y, 0) = sin(2πx)sin(2πy).

import numpy as np

h = .1
k = h/2
#h = .03125 # grid spacing on [0,1]^2, standard is h = .05
#k = (h**2)/4 # Time discretization
#k = h/5 
points = int(1/h) - 1 # Number of interior points along an axis
# Grid, u[0] -> [h,h], u[points^2 - 1] -> [1-h,1-h]
u = np.zeros(points**2)
# x and y direction discretization matrices for Laplacian
Dx = np.zeros((points**2,points**2))
Dy = np.zeros((points**2,points**2))

# Filling grid with initial conditions, i -> y, j -> x
for i in range(0,points):
    for j in range(0,points):
        u[i*points + j] = np.sin(2*np.pi*((j+1)*h))*np.sin(2*np.pi*((i+1)*h)) # ICs
        Dx[i*points + j,i*points + j] = -2/(h**2) # Diagonal 
        Dy[i*points + j,i*points + j] = -2/(h**2) # Diagonal

        if (j != 0):
            Dx[i*points + j, i*points + j - 1] = 1/(h**2) # Point to left
        if (j != points - 1):
            Dx[i*points + j, i*points + j + 1] = 1/(h**2) # Point to right
        if(i != 0):
            Dy[i*points + j, (i-1)*points + j] = 1/(h**2) # Point down
        if(i != points - 1):
            Dy[i*points + j, (i+1)*points + j] = 1/(h**2) # Point up

T = 1 # Final time
n = int(T/k) # number of time steps


DxInv = np.linalg.inv(np.identity(points**2) - (k/2)*Dx)
DyInv = np.linalg.inv(np.identity(points**2) - (k/2)*Dy)
Dx = np.identity(points**2) + (k/2)*Dx
Dy = np.identity(points**2) + (k/2)*Dy
M = np.linalg.multi_dot([DyInv,Dx,DxInv,Dy])

for i in range(0,n):
    #uStar = np.matmul(np.matmul(DxInv,Dy),u)
    #u = np.matmul(np.matmul(DyInv,Dx),uStar)
    u = np.matmul(M,u)
    
#for i in range(0,points):
#    for j in range(0,points):
#        print((j+1)*h, (i+1)*h, u[i*points + j])
        
error = 0
for i in range(0,points+2):
    for j in range(0,points+2):
        if ((i == 0 or j == 0) or (i == points + 1 or j == points + 1)):
            #print(j*h,i*h,0)
            continue
        else:
            error += np.abs(np.sin(2*np.pi*j*h)*np.sin(2*np.pi*i*h)*np.exp(-8*(np.pi**2)*T) - u[(i-1)*points + j-1])**2
            #print(np.sin(2*np.pi*j*h)*np.sin(2*np.pi*i*h)*np.exp(-8*(np.pi**2)*T), u[(i-1)*points + j-1])
            
            #print(j*h, i*h, u[(i-1)*points + j-1])

print(h,h*np.sqrt(error))
