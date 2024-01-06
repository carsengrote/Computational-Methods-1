import numpy as np

def v(x):
    return np.sin(np.pi*x/2)

def sinc(x):
    if (x == 0):
        return 1
    else:
        return np.sin(np.pi*x)/(np.pi*x)

def p(x,alpha):
    
    sum = 0
    sum = sum + v(0)*sinc(x)

    for m in range(1,alpha +1):

        sum = sum + v(m)*sinc(x - m)
        sum = sum + v(-m)*sinc(x + m)

    return sum

alphaValues = [1,2,4,8,16,32,64,128,256,512]

dx = .05 # For numerical integration, trapezoid
xValues = np.linspace(-5,5,int(10/dx))

#for x in xValues:
    #print(x,p(x,64),v(x))

for alpha in alphaValues:
    
    error = 0
    for x in xValues:
        error = error + ((p(x,alpha) - v(x))**2)*dx

    error = error - (dx/2)*((p(-5,alpha) - v(-5))**2 + (p(5,alpha)- v(5))**2)
    error = np.sqrt(error)
    print(alpha, error)
    

