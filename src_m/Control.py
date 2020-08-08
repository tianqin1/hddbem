# 8 11 15 16 19 23 26 27 29 31 36 43 106 107
# o  o  o  o  o o  o  o  o  x  x  x  x   x
import numpy as np
ng = 3
#ng =int(np.loadtxt(".ngg")[0])
# cplx=False
cplx = True
nregion = 3
latex = False
CBGXS = True
dataType = 'complex'
if(ng > 1 and cplx):
    cplx = True
    dataType = 'complex'
else:
    cplx = False
    dataType = 'float'
# general

nmed = nregion


width0 = 30
width1 = 15
width2 = 15
# width0 = 50
# width1 = 15
# width2 = 15
div0 = 500
div1 = 250
div2 = 250
# nmode = ng
# dl = width/nregion
# control parameters
omega = 1
eps = 1.0e-6
niter = 30
# # for medium control
# fuel = 0
# refl = 1

# flag control
pplot = False
converged = False
