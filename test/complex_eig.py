import numpy as np
aa = np.ones((3,3))
aa = np.array([[ 2 ,-3,  1], [ 1, -2,  1], [ 1 ,-3 , 2]])
bb,cc = np.linalg.eig(aa)
#print(bb)
#print(cc)


ini = np.zeros((3,3),dtype=complex)
ini2 = np.ones((3,3),dtype=complex)
#print(ini[0,0].real)
#print(ini,ini2)
c1  = np.complex(1,2)
c2 = np.complex(3,4)
div = c1/c2
print(c1)
print(c2)
print("c1-c2",c1-c2)
print(div*c2==c1)

nd = 2
cplx_mat_a = np.zeros((nd,nd),dtype=complex)
for i in range(nd):
	for j in range(nd):
		cplx_mat_a[i,j] = complex(np.random.random(1),np.random.random(1))

cplx_mat_b = np.zeros((nd,nd),dtype=complex)
for i in range(nd):
	for j in range(nd):
		cplx_mat_b[i,j] = complex(np.random.random(1),np.random.random(1))
print(cplx_mat_a*cplx_mat_b)
print(np.dot(cplx_mat_a,cplx_mat_b).real)
print(np.dot(cplx_mat_a.real,cplx_mat_b.real))