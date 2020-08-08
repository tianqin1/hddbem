import numpy as np
import math
import Control

# general green's function for 1d 
# k2:  buckling value (with square)
def gf(k2,x):
	x = abs(x)
	if Control.mg:
		k = np.sqrt(k2)
		return np.exp(1j*k*x)*1j/(2.0*k)
	if k2>0:
		k = np.sqrt(k2)
		return  -1/(2*k)*np.sin(k*x)
	else:# k2<0:
		k = np.sqrt(-k2)
#		if x == inf:
#			return 0
#		else:	
		return 1/(2*k)*np.exp(-k*x)

def bf(k2,x):
	x = abs(x)
	if Control.mg:
		k = np.sqrt(k2)
		return -np.exp(1j*k*x)/2.0
	if k2>0:
		k = np.sqrt(k2)
		return -0.5*np.cos(k*x)
	else:# k2<0:
		k = np.sqrt(-k2)
		return -0.5*np.exp(-k*x) 
			