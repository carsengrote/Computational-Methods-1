# Finite difference approximation for f''(x) given using x1 < x2 < x3 < x4

import numpy as np

def f(x):
	return np.exp(-1*x)*np.sin(x)

def f2(x):
	return -2*np.exp(-1*x)*np.cos(x)

def fpp(x1,x2,x3,x4):

	a = 2*(2*x2-x3-x4)/((x1-x2)*(x1-x3)*(x1-x4))
	b = 2*(2*x1-3*x2+x3+x4)/((x1-x2)*(x2-x3)*(x2-x4))
	c = 2*(x1-2*x2+x4)/((x1-x3)*(-1*x2+x3)*(x3-x4))
	d = 2*(x1-2*x2+x3)/((x1-x4)*(x2-x4)*(x3-x4))
	
	pp = a*f(x1) + b*f(x2) + c*f(x3) + d*f(x4)

	return pp

x1 = 0

for k in range(100,301):

	H = 10**(-1*k/100)
	x4 = H

	x2 = np.random.random()*H
	x3 = np.random.random()*H
	while(x3 == x2):
		x2 = np.random.random()*H
	
	if(x3 < x2):
		tmp = x2
		x2 = x3
		x3 = tmp
	
	guess = fpp(x1,x2,x3,x4)
	error = np.abs(guess-f2(x2))
	 
	print(H,error)
