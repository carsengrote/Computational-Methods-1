
import numpy as np
import numpy as np
import custom_plot as cplt

def timestep(u,uNext,k,h,b):
    
    for i in range(0,73):
        for j in range(0,160):

            if(vv[i,j] == 1): # Skip walls
                continue

            if ((36 <=i< 40) and (44 <=j< 48)): # Keep stinky pizza stinky
                uNext[i,j] = 1
                continue
            
            sum = 0
            if (i != 0): # UP
                if (vv[i-1,j] == 1):
                    sum = sum + u[i,j]
                else: 
                    sum = sum + u[i-1,j] 

            if(i != 72): # Down
                if(vv[i+1,j] == 1):
                    sum = sum + u[i,j]
                else:
                    sum = sum + u[i+1,j]
            
            if (j != 0): # Left
                if(vv[i,j-1] == 1):
                    sum = sum + u[i,j]
                else: 
                    sum = sum + u[i,j-1]

            if (j != 159): # Right
                if(vv[i,j+1] == 1):
                    sum = sum + u[i,j]
                else: 
                    sum = sum + u[i,j+1]

            sum = sum - 4*u[i,j]
            sum = sum*(k*b/(h**2))

            uNext[i,j] = sum + u[i,j]

    return uNext

vv = np.zeros((73,160)) # Map of the floor
file = open("van_vleck.txt", "r")
i = 0
for line in file:
    row = line.split()
    j = 0
    for num in row:
        vv[i,j] = int(num)
        j+=1
    i+=1

u = np.zeros((73,160)) # Numerical Solution, U_n
uNext = np.zeros((73,160)) # Solution at next timestep, U_{n+1}
for i in range(36,40):
    for j in range(44,48):
        u[i,j] = 1

h = 22.5 # cm, grid spacing
b = .55*(10000) # diffusion constant, .55 meters squared / second, converted to cm squared
k = (h**2)/(10*b) # time step
k = .01 # This is close, rounding down from given k
T = 100 # final time
n = int(T/k) # Number of time steps

#print("Time Step: ", k)
#print("Number of Steps: ", n)
#print("Final Time: ", n*k)

locations = np.matrix([[31,14],[58,103],[58,147]])
times = [0,0,0]
concentrations = np.zeros((n+1,3))

for i in range(1,n+1):

    u = timestep(u,uNext,k,h,b)

    for l in range(0,3):
        concentrations[n,l] = u[locations[l,0],locations[l,1]]
        if (u[locations[l,0],locations[l,1]]) >  10**(-4) and times[l] == 0:
            times[l] = i*k

    #print(i*k,concentrations[n,0],concentrations[n,1],concentrations[n,2])

u = u**(1/4)
print("C: ", times[0], " Q: ", times[1], " T: ", times[2])


# Printing out the image with Chris' code
# Load in the wall matrix
q=np.loadtxt("van_vleck.txt",dtype=np.int8)

# Call the first custom plotting routine
cplt.plot1("test1.png",u,q,0,1,4)

# Call the second custom plotting routine
#cplt.plot2("test2.png",u,q,-1.1,1.1,4)
