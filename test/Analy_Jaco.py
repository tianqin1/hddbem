import numpy as np
from Control import *

n=ng+1
jaco_1=np.loadtxt("jaco_1.txt",dtype=complex)
jaco_2=np.loadtxt("jaco_2.txt",dtype=complex)
diff=np.zeros((n*n),dtype=complex)
for i in range(n*n):
		diff[i]=jaco_1[i]-jaco_2[i]
		if(np.real(diff[i])<1.e-7):
			diff[i]=complex(0,np.imag(diff[i]))
		if(np.imag(diff[i])<1.e-7):
			diff[i]=complex(np.real(diff[i]),0)

#print(jaco_1.reshape((n,n)))			
print(diff.reshape((n,n)))