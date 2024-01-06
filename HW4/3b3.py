import numpy as np
import sys

mValues = [200,400,800,1600,3200,6400,12800] # Number of domains
for mb in mValues:
    m = int(mb)
    h = 2*np.pi/m # Space discretization
    c = 10/3 # CFL condition
    k = h/(3*c) # Time discretization
    T = 3*np.pi/np.sqrt(5) # Longest time

    def A(i,h):
        return 2 + (4/3)*np.sin((i+1/2)*h)

    def Qnext(Q,Ai,Ah,h,k):
        Qcenter = Q
        Qleft = np.roll(Q,1)
        Qright = np.roll(Q,-1)
        Acenter = Ai
        Aleft = np.roll(Ai,1)
        Aright = np.roll(Ai,-1)
        ArightHalf = np.roll(Ah,-1)

        Fleft  = (Aleft*Qleft + Acenter*Qcenter)/2 - (Ah*k/h)*(Acenter*Qcenter - Aleft*Qleft)/2
        Fright = (Acenter*Qcenter + Aright*Qright)/2 - (ArightHalf*k/h)*(Aright*Qright - Acenter*Qcenter)/2

        return Qcenter-(k/h)*(Fright - Fleft)

    finalT = T
    n = int(finalT/k) # number of steps
    k =finalT/n # Time discretization

    q = np.zeros(m) # Initial concentration on intervals
    qNext = np.zeros(m)
    Ah = np.zeros(m)
    Ai = np.zeros(m)
    for i in range(0,m):
        q[i] = np.exp(np.sin((i+1/2)*h)+(1/2)*np.sin(4*(i+1/2)*h))
        #q[i] = np.max(np.pi/2 - np.abs((i + 1/2)*h - np.pi),0)
        Ah[i] = A(i-1/2,h)
        Ai[i] = A(i,h)

    for j in range(0,n):
        q = Qnext(q,Ai,Ah,h,k)

    sum = 0
    for i in range(0,m):
        #print((i+1/2)*h, q[i])
        sum = sum + h*(((np.exp(np.sin((i+1/2)*h)+(1/2)*np.sin(4*(i+1/2)*h))) - q[i])**2)
	#sum = sum + h*((np.max(np.pi/2 - np.abs((i+1/2)*h - np.pi),0) - q[i])**2)

    error = np.sqrt(sum)
    print(m,error)
