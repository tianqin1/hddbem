import numpy as np
import matplotlib.pyplot as plt
from Material import Material
import Green
from Control import *
import Control as ctrol


if(ng>3 and cplx):
	fuel = Material()
	file1 = "../xs/cplx_xs/"+str(ng)+"g/cplx_1d_xs.txt"
	#print(str(file1))
	readxs1d = np.reshape(np.loadtxt(file1),(2,ng,6))
	# 2 region, 3 group, 4 kind of 1d xs 
	xs1d=np.zeros((2,ng,4))
	xs1d[:,:,0]=readxs1d[:,:,3]
	xs1d[:,:,1]=readxs1d[:,:,2]
	xs1d[:,:,2]=readxs1d[:,:,1]
	xs1d[:,:,3]=readxs1d[:,:,-1]
	file2="../xs/cplx_xs/"+str(ng)+"g/cplx_2d_xs.txt"

	xs2d = np.reshape(np.loadtxt(file2),(2,ng,ng))

	fuel.set(xs1d[0])
	fuel.set_smat(xs2d[0])
	refl = Material()
	refl.set(xs1d[1])
	refl.set_smat(xs2d[1])

# ### 1g xs
elif(ng==1 and not cplx):
	fuel = Material()
	fuel.set([[1.36, 0.0181, 0.0279, 1.0]])
	#fuel.set_smat( [[0.0, 0.0178], [0.0, 0.0]])
	refl = Material()
	#refl.set([[0.55, 0.0181, 0.0, 1.0]])
	refl.set([[0.55, 0.0127, 0.0, 1.0]])
	#refl.set_smat( [[0.0, 0.0476], [0.0, 0.0]])

elif(ng==2 and not cplx):
### 2g xs
	fuel = Material()
	fuel.set([[1.58, 0.0032, 0.0, 1.0],[0.271, 0.0930, 0.168, 0.0]])
	fuel.set_smat( [[0.0, 0.0178], [0.0, 0.0]])
	refl = Material()
	refl.set([[1.41, 0.0, 0.0, 1.0],[0.117, 0.0191, 0.0, 0.0]])
	refl.set_smat( [[0.0, 0.0476], [0.0, 0.0]])
elif(ng==3 and not cplx):
# ### 3g xs real
	fuel = Material()
	fuel_xs1d=[[1.60088553e+00, 4.14769275e-03, 4.63179141e-03, 9.99989828e-01],
			  [5.40995656e-01, 3.67921129e-02, 2.45307219e-02, 5.41253125e-06],
			  [3.19523666e-01,1.16339868e-01,1.94546876e-01,1.23390469e-10]]

	fuel.set(fuel_xs1d)
	fuel.set_smat( [[4.54551305e-01, 3.28783130e-02, 2.57438338e-05],
	                [0.00000000e+00, 8.11562493e-01, 7.85075787e-02],
	                [0.00000000e+00, 3.67013825e-03, 1.42044709e+00]])

	refl = Material()
	refl.set([[1.94496989e+00, 2.10489436e-04, 0.0, 3.64494000e-01],
			  [5.74634121e-01, 1.48218560e-03, 0.0, 3.17764000e-01],
	          [4.45679857e-01,1.87048450e-02, 0.0,3.17764000e-01 ]])
	refl.set_smat( [[6.06991999e-01, 7.10768286e-02, 5.63638376e-05],
	                [0.00000000e+00, 1.32556054e+00, 1.89770969e-01],
	                [0.00000000e+00, 1.90967691e-04, 3.22053777e+00]])

elif(ng==3 and cplx):
	fuel = Material()
	fuel_xs1d=[[1.40941976e+00, 1.07180561e-02, 7.56465032e-03, 9.99995238e-01],
     [6.40003843e-01, 2.29853993e-02, 3.16497075e-02, 2.86925000e-09],
     [3.21927345e-01, 1.15174972e-01, 1.92555606e-01, 1.23390469e-10]]

	fuel.set(fuel_xs1d)
	fuel.set_smat( [[5.37345373e-01, 1.77258266e-02, 2.18569668e-03],
     [8.29957684e-05, 6.48012357e-01, 2.44437484e-01],
     [0.00000000e+00, 3.82527695e-03, 1.41148228e+00]])

	refl = Material()
	refl.set([[1.54025533e+00, 2.24020966e-04, 0.00000000e+00, 5.51414000e-01],
     [5.69188583e-01, 3.21135432e-03, 0.00000000e+00, 1.30844000e-01],
     [4.71488040e-01, 1.57704709e-02, 0.00000000e+00, 3.17764000e-01]])
	refl.set_smat([[8.04398991e-01, 3.79949212e-02, 5.08707986e-03],
     [0.00000000e+00, 1.04632637e+00, 5.16367152e-01],
     [0.00000000e+00, 8.94450256e-04, 2.90370915e+00]])


med = np.array([fuel,refl])
print(fuel.show_xs())
print(fuel.show_scatter_xs())
print(refl.show_xs())
print(refl.show_scatter_xs())


def bem(xinp):
	if(cplx):
		k_ini=xinp[0].real
	else:
		k_ini=xinp[0]

	inmode=xinp[1:]
	mode0 = np.zeros(ng,dtype=dataType)
	inmode_r = np.zeros(ng,dtype=dataType)
	dmodex = np.zeros(ng,dtype=dataType)
	dmodex_r = np.zeros(ng,dtype=dataType)
	amat = np.zeros((nmed,ng,ng))
	bmat = np.zeros((nmed,ng),dtype=dataType)
	cmat = np.zeros((nmed,ng,ng),dtype=dataType)

	for i in range(nmed):
		amat[i] = med[i].calc_matrix_a(k_ini)
		bmat[i], cmat[i] = np.linalg.eig(amat[i])
	for i in range(nmed):
		for j in range(ng-1):
			cmat[i,j,:]=cmat[i,j,:]/cmat[i,-1,:]
		cmat[i,-1,:]=np.ones(ng)

	dmat = np.zeros((nmed,ng,ng))
	for i in range(nmed):
		for g in range(ng):
			dmat[i,g,g] = med[i].dif(g)

	fuel = 0
	refl = 1
	# neutron flux of fuel reguon at inner boundary
	inflx=np.zeros(ng,dtype=dataType)
	#if(ng>1):
	inmode_r = np.dot(np.linalg.inv(cmat[1]),np.dot(cmat[0],inmode))
	#else:
	#inmode_r=(cmat[0]*inmode)/cmat[1]
	for j in range(ng):
		b2 = bmat[fuel][j]
		tmp1 = inmode[j]
		tmp2=0
		tmp=0
		if (b2>0):
			#print('b2>0 ? b2=',b2)
			b=np.sqrt(b2)
			#print('b2>0? b=',b)
			tmp2 = (tmp1/np.cos(b*width))
			tmp=(-tmp2+np.cos(b*width)*tmp1)*b/(np.sin(b*width))
		else:
			b=np.sqrt(-b2)
			#print("### b = ",b)
			tmp2=(tmp1*2.0/(np.exp(b*width)+np.exp(-b*width)))
			tmp=b*tmp1-tmp2*np.exp(-b*width)
		#print(tmp2)
		mode0[j]=tmp2
		dmodex[j]=tmp
# for reflector region with infinite length 
	for j in range(ng):
		k2=bmat[refl]
		k=np.sqrt(-k2)
		tmp=-k*inmode_r[j]
		dmodex_r[j]=tmp[j]

	curr_r = np.zeros(ng,dtype=dataType)
	curr_f = np.zeros(ng,dtype=dataType)
	#for i in range(2):
	#print(dmodex)
	#print(dmodex_r)
	# if(ng>1):
	curr_f = np.dot(dmat[fuel],np.dot(cmat[fuel],dmodex))
	curr_r = np.dot(dmat[refl],np.dot(cmat[refl],dmodex_r))
	# else:
	# 	curr_r = dmodex_r
	# 	curr_f = dmodex

	curr_f = curr_f -curr_r
	#print("curr_f",curr_f)
	oup=np.zeros(ng+1,dtype=dataType)
	norm = 1-np.sum(inmode)
	for i in range(ng+1):
		if(i==0):
			oup[i]=norm
		else:
			oup[i]=curr_f[i-1]
#	oup = np.array([norm,curr_f[0],curr_f[1],curr_f[2]])
	#print(Control.converged)
	if(ctrol.converged):
		print("==>> Saving converged data ...")
		np.savetxt('.mode0.txt',mode0)
		np.savetxt('.inmode.txt',inmode)
		np.savetxt('.dmodex.txt',dmodex)
		np.savetxt('.inmode_r.txt',inmode_r)
		np.savetxt('.dmodex_r.txt',dmodex_r)
		np.savetxt('.amat.txt',np.reshape(amat,(nmed*ng*ng)))
		np.savetxt('.bmat.txt',bmat)
		np.savetxt('.cmat.txt',np.reshape(cmat,(nmed*ng*ng)))
	return oup

def get_flx(inp):
	# if(not pplot):
	# 	print("Please turn on the pplot flag on ...")
	# else:

	bem(inp)
	print("==>> Reconstrcuting the neutron flux ...")

	mode0=np.loadtxt('.mode0.txt',dtype=dataType)
	inmode = np.loadtxt('.inmode.txt',dtype=dataType)

	dmodex = np.loadtxt('.dmodex.txt',dtype=dataType)
	inmode_r = np.loadtxt('.inmode_r.txt',dtype=dataType)
	dmodex_r = np.loadtxt('.dmodex_r.txt',dtype=dataType)
	bmat = np.loadtxt('.bmat.txt',dtype=dataType)
	cmat = np.reshape(np.loadtxt('.cmat.txt',dtype=dataType),(nmed,ng,ng))
	div = 500
	dx = width/div
	x= dx*0.5
	#cmat_f = cmat[0]
	#print(bmat)
	#print(cmat[0])
	mode_tmp_f = np.zeros(ng,dtype=dataType)
	tmp=0
	xflx_f = np.zeros((div,ng+1),dtype=dataType)
	xmodeflx_f = np.zeros((div,ng+1),dtype=dataType)
### one group neutron flux
	if(ng==1):
		# fuel region
		for i in range(div):
			for j in range(ng):
				b2=bmat[0]
				if(b2>0):
					b=np.sqrt(b2)
					tmp=-0.5/b*np.sin(b*(width-x))*dmodex
					tmp+=0.5*np.cos(b*(width-x))*inmode+0.5*np.cos(-b*x)*mode0
					mode_tmp_f =tmp
				else:
					b=np.sqrt(-b2)
					tmp=0.5/b*np.exp(-b*(width-x))*dmodex
					tmp+=0.5*np.exp(-b*(width-x))*inmode+0.5*np.exp(-b*x)*mode0
					mode_tmp_f=tmp
				print("mode_tmp_f",mode_tmp_f)
			xmodeflx_f[i,0]=x
			xmodeflx_f[i,1:ng+1] = mode_tmp_f
			xflx_f[i,0]=x
			xflx_f[i,1:ng+1] = np.dot(cmat[0],mode_tmp_f)
			x+=dx
		#print(xmodeflx_f)
		x=dx*0.5
		tmp=0
		xflx_r = np.zeros((div,ng+1),dtype=dataType)
		xmodeflx_r = np.zeros((div,ng+1),dtype=dataType)
		mode_tmp_r = np.zeros(ng,dtype=dataType)
		# reflector retion
		for i in range(div):
			for j in range(ng):
				b2=bmat[1]
				b = np.sqrt(-b2)
				tmp = -0.5/b*np.exp(-b*x)*dmodex_r
				tmp+=0.5*np.exp(-b*x)*inmode_r
				mode_tmp_r = tmp
			xmodeflx_r[i,0]=x+width
			xmodeflx_r[i,1:ng+1] = mode_tmp_r
			xflx_r[i,0]=x+width
			xflx_r[i,1:ng+1]=np.dot(cmat[1],mode_tmp_r)
			x+=dx
	else:
###=============== multi group ================
		# fuel region
		for i in range(div):
			for j in range(ng):
				b2 = bmat[0,j]
				#print("b2",b2)
				if(b2>0):
					#print("b2>0",b2)
					b=np.sqrt(b2)
					#print("b2>0, b=",b)
					#tmp=-0.5/b*np.sin(b*(width-x))*dmodex[j]
					#tmp+=0.5*np.cos(b*(width-x))*inmode[j]+0.5*np.cos(-b*x)*mode0[j]
					#tmp=
					tmp=-0.5/b*np.sin(b*(width-x))*dmodex[j]+0.5*np.cos(b*(width-x))*inmode[j]+0.5*np.cos(-b*x)*mode0[j]
					mode_tmp_f[j] =tmp
				else:
					#print("b2 ",b2)
					b=np.sqrt(-b2)
					#tmp=0.5/b*np.exp(-b*(width-x))*dmodex[j]
					tmp=0.5/b*np.exp(-b*(width-x))*dmodex[j]+0.5*np.exp(-b*(width-x))*inmode[j]+0.5*np.exp(-b*x)*mode0[j]
					mode_tmp_f[j]=tmp

			xmodeflx_f[i,0]=x
			xmodeflx_f[i,1:ng+1] = mode_tmp_f
			xflx_f[i,0]=x
			xflx_f[i,1:ng+1] = np.dot(cmat[0],mode_tmp_f)
			x+=dx
		#print(xmodeflx_f)
		x=dx*0.5
		tmp=0
		xflx_r = np.zeros((div,ng+1),dtype=dataType)
		xmodeflx_r = np.zeros((div,ng+1),dtype=dataType)
		mode_tmp_r = np.zeros(ng,dtype=dataType)
		# reflector retion
		for i in range(div):
			for j in range(ng):
				if(ng==1):
					b2=bmat[1]
				else:
					b2 = bmat[1,j]
				b = np.sqrt(-b2)
				tmp = -0.5/b*np.exp(-b*x)*dmodex_r[j]
				tmp+=0.5*np.exp(-b*x)*inmode_r[j]
				mode_tmp_r[j] = tmp
			xmodeflx_r[i,0]=x+width
			xmodeflx_r[i,1:ng+1] = mode_tmp_r
			xflx_r[i,0]=x+width
			xflx_r[i,1:ng+1]=np.dot(cmat[1],mode_tmp_r)
			x+=dx
	#print(xmodeflx_r)
	flx=np.array([xflx_f,xflx_r])
	mflx = np.array([xmodeflx_f,xmodeflx_r])
	#print(flx.shape)
	#pplot=True
	plt.xlabel("Distance from the core center ${x} [m]$")

	flux = np.concatenate((flx[0],flx[1]),axis=0)
	mflux = np.concatenate((mflx[0],mflx[1]),axis=0)
	print("==>> Saving the neutron flux and mode neutron flux ...")
	np.savetxt(".flux.txt",flux)
	np.savetxt(".mode.flux.txt",mflux)
	#print(mflx.shape)
	print("==>> Plotting the neutron flux, [Figure 1]... ")
	plt.figure(1)
	plt.title("Neutron flux distribution of 1-D half core")

	for i in range(ng):
		plt.plot(flux[:,0].real,flux[:,i+1].real,label="Group " +str(i+1))
	plt.xlabel("Distance from the core center ${x} [m]$")

	plt.grid( which='both',axis="both")
	plt.ylabel("Neutron Flux ${\phi}$[x]")


	print("==>> Plotting the mode neutron flux, [Figure 2] ...")
	plt.figure(2)

	for g in range(ng):
		plt.subplot(ng,nmed,g*nmed+1)
		if(g==0):
			plt.title("Core region")
		for mode in range(ng):
			if(mode==1):
				plt.plot(mflx[0,:,0].real,np.real(cmat[0,g,mode]*mflx[0,:,mode+1]),".-",linewidth=0.7,markersize=1,label="Mode "+str(mode+1))
			else:
				plt.plot(mflx[0,:,0].real,np.real(cmat[0,g,mode]*mflx[0,:,mode+1]),".-",linewidth=0.7,markersize=1,label="Mode "+str(mode+1))

		plt.plot(mflx[0,:,0].real,np.real(np.dot(np.transpose(cmat[0,g,:]),np.transpose(mflx[0,:,1:ng+1]))),"k-",markersize=1,label="Combination")
		plt.grid( which='both',axis="both")
		if(ng<5):
			plt.legend(loc='lower left',prop={'size': 7})
		plt.xlabel(" ${x} [m]$")
		plt.ylabel("Group "+str(g+1))

		plt.subplot(ng,nmed,g*nmed+2)
		if(g==0):
			plt.title("Reflector region")
		for mode in range(ng):
			plt.plot(mflx[1,:,0].real,np.real(cmat[1,g,mode]*mflx[1,:,mode+1]),".-",linewidth=0.7,markersize=1,label="Mode "+str(mode+1))
		plt.plot(mflx[1,:,0].real,np.real(np.dot(np.transpose(cmat[1,g,:]),np.transpose(mflx[1,:,1:ng+1]))),"k-",markersize=1,label="Combination")
		plt.xlabel(" ${x} [m]$")
		plt.grid( which='both',axis="both")
		if(ng<5):
			plt.legend(loc='upper right',prop={'size': 7})
	plt.show()
	print("===================================================")
	if(ng<15):
		print("==>> Printing buckling value")
		for i in range(2):
			if(i==0):
				print("Core ",end="")
			else:
				print("Reflector ",end="")
			for j in range(ng):
				if(cplx):
					print("&", "{0.real:7.5e}{0.imag:+7.5e}j".format(bmat[i,j]),end="")
				else:
					print("&", "%7.5e"%bmat[i,j],end="")
			print("\\\\")
		print("==>> Printing cmat value")	
		for i in range(ng):
			for j in range(ng):
				if(cplx):
					print("&", "{0.real:7.5e}{0.imag:+7.5e}j".format(cmat[0,i,j]),end="")
				else:
					print("&", "%7.5e"%cmat[0,i,j],end="")
			print("\\\\")
		for i in range(ng):
			for j in range(ng):
				if(cplx):
					print("&", "{0.real:7.5e}{0.imag:+7.5e}j".format(cmat[1,i,j]),end="")
				else:
					print("&", "%7.5e"%cmat[1,i,j],end="")
			print("\\\\")
	return flux

def get_inmode(inp_ngg,inp_keff):
	ngg=inp_ngg
	keff=inp_keff 
	nmesh=1000
	k_ini=keff
	#flux=np.loadtxt('./cplx_xs/'+str(ngg)+'g/flux_cbz.txt',dtype=dataType)
	flux=np.loadtxt('./cplx_xs/'+str(ngg)+'g/flux.txt',dtype=dataType)
	#print(flux.shape)
	#print(flux[:,0].shape)
	flx_tmp=flux[:,1:]
	flx=np.zeros((nmesh,ngg),dtype=dataType)
	flx=flx_tmp
	#for g in range(ngg):
	#	for i in range(nmesh):
	#		flx[i,g]=flx_tmp[nmesh*g+i]
	
	flx0=flx[0,:]
	flx_x=flx[499,:]
	flx_xr=flx[500,:]
	#flx_x=(flx_x+flx_xr)/2.0
	#print(flx_x)

	inmode=np.zeros(ngg,dtype=dataType)
	mode0 = np.zeros(ngg,dtype=dataType)
	inmode_r = np.zeros(ngg,dtype=dataType)
	dmodex = np.zeros(ngg,dtype=dataType)
	dmodex_r = np.zeros(ngg,dtype=dataType)
	
	amat = np.zeros((nmed,ngg,ngg))
	bmat = np.zeros((nmed,ngg),dtype=dataType)
	cmat = np.zeros((nmed,ngg,ngg),dtype=dataType)

	for i in range(nmed):
		amat[i] = med[i].calc_matrix_a(k_ini)
		bmat[i], cmat[i] = np.linalg.eig(amat[i])
	for i in range(nmed):
		for j in range(ngg-1):
			cmat[i,j,:]=cmat[i,j,:]/cmat[i,-1,:]
		cmat[i,-1,:]=np.ones(ngg)

	dmat = np.zeros((nmed,ngg,ngg))
	for i in range(nmed):
		for g in range(ngg):
			dmat[i,g,g] = med[i].dif(g)
	#np.dot(inmode,cmat[0])=flx
	#print(cmat[0].shape)
	inmode=np.dot(np.linalg.inv(cmat[0]),flx_x)
	#print(inmode)
	return inmode

## for debug
if __name__ == '__main__':
	if(ng==2):
		keff=1.3588487076917326
		#keff=1.3
		inmode = np.ones(ng,dtype=dataType)
		inp = np.ones(ng+1)
		inmode=np.array([0.3592456,0.6407544,0.5],dtype=dataType)
		for i in range(ng+1):
			if(i==0):
				inp[i]=keff
			else:
				inp[i]=inmode[i-1]
		print(inp)
		#print(bem(inp))
		inp = bem(inp)
		print(inp)
	elif(ng==1):
		inikeff=1.4110198630169637
		inmode=1
		inp=np.zeros(ng+1)
		inp[0]=inikeff
		inp[1]=inmode
		print(bem(inp))
	elif(ng==26):
		inmode=np.array([    (4.092329189924370225e-01-2.077386268162504223e-14j),
							 (9.589616304541248581e-02+1.865482333924856710e-13j),
							 (8.259992446648199704e-02+2.453077998118275522e-13j),
							 (2.842810998476946738e-01+1.619100821271635566e-13j),
							 (1.292258688685308209e-02+2.236711413825824792e-13j),
							 (6.404222749752655139e-03-2.159632179865490428e-02j),
							 (6.404222750291950310e-03+2.159632179863119061e-02j),
							 (-3.905469005985604761e-03-1.265853678475723464e-02j),
							 (-3.905469005664603289e-03+1.265853678481873058e-02j),
							 (1.756804355688295488e-02-2.416560310780374934e-02j),
							 (1.756804355744772186e-02+2.416560310763356950e-02j),
							 (-2.364299583985689929e-02+1.177280080071516598e-02j),
							 (-2.364299584009146513e-02-1.177280080045768618e-02j),
							 (4.854237309868530102e-03+1.382107331579375976e-03j),
							 (4.854237309907380969e-03-1.382107331663969549e-03j),
							 (2.586729911886907415e-02-1.603217671763968941e-02j),
							 (2.586729911920973568e-02+1.603217671722353271e-02j),
							 (6.324282170772228817e-02-4.209875507433057968e-13j),
							 (-2.019408896812561105e-04+1.652809757433005765e-15j),
							 (7.263424983280934635e-04-1.003722844039829379e-14j),
							 (-5.505051917566706164e-04+3.627760181362434148e-15j),
							 (1.177575383623623848e-04-1.398539233794634127e-15j),
							 (5.528018403427196808e-04+8.466622525333136404e-03j),
							 (5.528018402109014049e-04-8.466622525349624950e-03j),
							 (-1.010850026066778930e-03+5.882032591298134693e-15j),
							 (-2.652598336972353168e-03+1.664569845068695286e-14j)])
		inp=np.zeros(ng+1,dtype=dataType)
		inp[0]=complex(1.22268727,0)
		inp[1:]=inmode
		print(bem(inp))
	elif(ng==29):
		inmode=np.array([ (3.646004205262999287e-01-8.163389871816788013e-17j),
							 (8.587852548644558937e-02-1.197938074220416344e-16j),
							 (7.705296812203533963e-02-1.410732712723418675e-16j),
							 (3.185421038426414220e-02-7.508947982865181696e-17j),
							 (5.388877274765645903e-02-1.433244252522904594e-16j),
							 (3.506410100447369255e-02-9.951673571528886191e-17j),
							 (2.522702828516211124e-02-7.766575520950835441e-17j),
							 (1.152537475222926988e-02-3.438109441506454807e-17j),
							 (1.791161683427161533e-02-5.333770322401675029e-17j),
							 (6.764083245745410033e-02+3.790314631825188076e-16j),
							 (-1.653477090265554478e-02-9.842501952066643340e-17j),
							 (2.978312400587038167e-02-5.327538144901370135e-03j),
							 (2.978312400587021513e-02+5.327538144901924379e-03j),
							 (2.521337178736517137e-02-1.768719884499605479e-04j),
							 (2.521337178736463708e-02+1.768719884503364408e-04j),
							 (2.276177644402364864e-02+1.878772115000926918e-02j),
							 (2.276177644402358966e-02-1.878772115000935938e-02j),
							 (-3.383577845950153096e-02-2.401851056642379983e-16j),
							 (3.862583020149108759e-02+3.973409639923454535e-16j),
							 (-2.144694598778851095e-03-6.697269196458619007e-17j),
							 (3.179657088104104537e-03-1.731216894588276457e-16j),
							 (9.389015786373954356e-03-2.840902670356344630e-17j),
							 (6.831457749685936574e-03-2.040318924880987894e-17j),
							 (1.113123491062486571e-02+6.859421029336063458e-18j),
							 (7.830649281885765020e-03-3.463920534384341781e-17j),
							 (1.341474408717098499e-02-3.697068559633246145e-17j),
							 (1.541486360439146483e-02-4.359799297222251633e-17j),
							 (9.092491203999426147e-03-2.317133841835882798e-17j),
							 (1.144490497239842274e-02-2.996195270482481645e-17j)
							])
		inp=np.zeros(ng+1,dtype=dataType)
		inp[0]=complex(1.21638242,0)
		inp[1:]=inmode
		print(bem(inp))

	elif(ng==107):
		keff=1.21763 
		inmode=get_inmode(107,keff)
		inp=np.ones(ng+1,dtype=dataType)
		for i in range(1,ng+1):
			inp[i]=inmode[i-1]

		inp[0]=complex(keff,0)
		print(bem(inp))
	else:
		print("Debug", ng,"Group ")
		keff=1.216923684972949
		inmode = np.ones(ng,dtype=dataType)/ng
		if(cplx):
			for i in range(ng):
				inmode[i]=complex(1,1)
		inp = np.ones(ng+1,dtype=dataType)
		inp[0]=complex(keff,1.0)
		#get_inmode()
		#bem(inp)
	inp_ngg=29
	inp_keff=1.21638242
	print(get_inmode(inp_ngg,inp_keff))
