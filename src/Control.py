# 8 11 15 16 19 23 26 27 29 31 36 43 106 107
# o  o  o  o  o o  o  x  o  x  x  x  x   x
import numpy as np
#ng =int(np.loadtxt("ngg")[0])
ng=3
#cplx=False
#sp3=False
cplx=True
dataType='complex'
if(ng>2 and cplx):
	cplx=True
	dataType='complex'
else:
	cplx=False
	dataType='float'
# general
nregion = 2
nmed = 2
nnode = 4
width = 30
nmode = ng
dl = width/nregion
# control parameters
alpha = 1.0
eps = 1.0e-8
niter = 30
# # for medium control
# fuel = 0
# refl = 1

# flag control
pplot=False
converged = False
