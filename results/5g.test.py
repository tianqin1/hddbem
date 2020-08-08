import numpy as np
nmed=2
ng=5
dataType=complex
cmat=np.reshape(np.loadtxt('5g.cmat.txt',dtype=dataType),(nmed,ng,ng))
bmat=np.loadtxt('5g.bmat.txt',dtype=dataType)
inmode=np.loadtxt('5g.inmode.txt',dtype=dataType)
flux=np.loadtxt('5g.flux.txt',dtype=dataType)

print(np.dot(cmat[0],inmode))