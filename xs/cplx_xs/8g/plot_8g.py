import numpy as np
import matplotlib.pyplot as plt
import os
ng = 8
cmd = 'cp /Users/tian/Documents/research/bem/.flux.txt /Users/tian/Documents/research/bem/cplx_xs/'+str(ng)+'g'+'/flux.txt'
os.system(cmd)
flux = np.loadtxt("flux_cbz.txt")
bem_flux=np.loadtxt("flux.txt",dtype=complex)
nmesh=1000
flx=np.zeros((nmesh,ng))
x= np.zeros(nmesh)
x=flux[0:nmesh,0]

for i in range(nmesh):
    for g in range(ng):
        flx[i,g]=flux[nmesh*g+i,1]


a=np.zeros(ng)
for i in range(ng):
	a[i] = flx[0,i]/bem_flux[0,i+1].real
	plt.plot(x,flx[:,i]/a[i],'k-',linewidth=1.3,label='FVM'+str(i+1)+'G')

colors=['red','magenta','purple','blue','green','blue','blue','blue']
for x in range(0,nmesh,50):
	for i in range(ng):
		plt.plot(bem_flux[x,0].real,bem_flux[x,i+1].real,'-s',linewidth=1.3,color=colors[i])
#		plt.plot(bem_flux[x,0].real,bem_flux[x,i+1].real,'-s',linewidth=1.3)

#plt.legend(['FVM 1G','FVM 2G','FVM 3G','FVM 4G','FVM 5G','HDD-BEM 1G','HDD-BEM 2G','HDD-BEM 3G','HDD-BEM 4G','HDD-BEM 5G'])
plt.grid()
plt.xlabel("Distance from core center x[m]")
plt.ylabel("Neutron flux distribution ${\phi}[/(m^2*s)$]")
plt.show()
