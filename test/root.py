import cmath
import numpy as np

import random
from scipy import optimize

cplx=True

if(cplx):
	dataType = 'complex'
else:
	dataType = 'float'
size=3

# def function_fx(x1,x2,x3):
# 	x=np.zeros(size,dtype=dataType)
# 	x[0]=x1
# 	x[1]=x2
# 	x[2]=x3
# 	# size = np.shape(x)[0]
# 	y = np.zeros(size, dtype=dataType)
# 	if(cplx):
# 	    y[0] = x[0] * x[1] * x[2] - complex(6,0);
# 	    y[1] = x[0] * x[0] * x[1] +  x[1] * x[1] * x[2] +  x[2] * x[2] * x[0] - complex(23,0);
# 	    y[2] = np.exp(x[0] + x[1] + x[2]) - complex(403,0);
# 	else:
# 	    y[0] = x[0] * x[1] * x[2] - 6;
# 	    y[1] = x[0] * x[0] * x[1] + x[1] * x[1] * x[2] + x[2] * x[2] * x[0] - 23;
# 	    y[2] = np.exp(x[0] + x[1]) - x[2];
# 	return y[0],x[1],y[2]

def function_fx(x):
	size = np.shape(x)[0]
	y = np.zeros(size, dtype=dataType)
	if(cplx):
	    y[0] = x[0] * x[1] * x[2] - complex(6,0);
	    y[1] = x[0] * x[0] * x[1] +  x[1] * x[1] * x[2] +  x[2] * x[2] * x[0] - complex(23,0);
	    y[2] = np.exp(x[0] + x[1] + x[2]) - complex(403,0);
	else:
	    y[0] = x[0] * x[1] * x[2] - 6;
	    y[1] = x[0] * x[0] * x[1] + x[1] * x[1] * x[2] + x[2] * x[2] * x[0] - 23;
	    y[2] = np.exp(x[0] + x[1] + x[2]) - 403;
	return y
np.random.seed(1)
x0=np.random.random(3)
print(x0)
# method = 
# hybr				not working
# lm                not working
# broyden1			not working
# broyden2			sometimes works but not converge
# anderson			not working
# linearmixing		works but never converge
# diagbroyden 		not working
# excitingmixing 	works but never converge
# krylov 			not working
# df-sane 	 	 	works and sometimes converge
my_options={"line_search":"cruz"}
sol = optimize.root(function_fx, x0,method='df-sane',options=my_options)
print(sol)
#1.30582274, 1.39181873, 3.30129509
#0.99953451, 2.0040086 , 2.99539345
