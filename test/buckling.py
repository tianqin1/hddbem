import numpy as np
bmat=np.loadtxt("../src/.bmat.txt",dtype=complex)
cmat=np.loadtxt("../src/.cmat.txt",dtype=complex)
cmat=np.reshape(cmat,(bmat.shape[0],bmat.shape[1],bmat.shape[1]))
# for bare slab reactor of fuel region
B2 = (np.pi/30)
print(B2)

cb=np.zeros((bmat.shape[0],bmat.shape[1]),dtype=complex)
for i in range(bmat.shape[0]):
	cb[i]=np.dot(cmat[i],bmat[i])
print(sum(cb[0]))