import numpy as np
from Control import *
from scipy.linalg import ldl

amat = np.reshape(np.loadtxt('.amat.txt',dtype=dataType),(nmed,ng,ng))
a=amat[0]
#print(a)
lu, d, perm = ldl(a, lower=0)
print(lu)
print(d)
print(perm)