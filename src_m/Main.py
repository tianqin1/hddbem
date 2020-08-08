import cmath
import numpy as np
# import numdifftools as nd
from scipy import optimize
from mpmath import findroot
import os
import Control as ctrol
from Control import *
import HDDBEM as hddbem
#from TwoGroup import *


def hessian ( x):
    N = x.size
    h = np.zeros((N,N),dtype=dataType)
    df_0 = hddbem.bem ( x )[1]
    for i in range(N):
        xx0 = 1.*x[i]
        x[i] = xx0 + complex(eps,eps)
        df_1 = hddbem.bem ( x )[1]
        h[i,:] = (df_1 - df_0)/complex(eps,eps)
        x[i] = xx0
    return h
def fdjac(x):
#calculate jacobian matrix for real value function
	n = np.shape(x)[0]
	xh=np.zeros(n)
	for i in range(n):
		xh[i] = x[i]
	df = np.ones((n,n))
	fvec = np.zeros(n)
	#print(x)
	fvec = hddbem.bem(x)
	#print(fvec)
	temp = 0
	h = 0
	for j in range(n):
		#print("mode: ",j)
		temp = xh[j]
		h = eps * abs(temp)
		if abs(h)==0.0:
			h = eps
		xh[j] = temp + h
		h = xh[j] - temp
		f = np.zeros(n)
		#print(xh)
		f = hddbem.bem(xh)
		#print(f)
		for i in range(n):
			#print(((f[i]-fvec[i])/h))
			df[i][j] = (f[i]-fvec[i])/h
		xh[j] = temp

	return df

def fdjac2(x):
#calculate jacobian matrix for complex value function
	n = np.shape(x)[0];
	xh = np.zeros(n,dtype=dataType)
	for i in range(n):
		xh[i] = x[i]
	df = np.ones((n,n), dtype=dataType);
	for i in range(n):
		for j in range(n):
			df[i,j]=np.complex(1,1)
	#print(df)
	fvec = np.zeros(n, dtype=dataType)
	#print(x)
	fvec = hddbem.bem(x)
	temp = np.complex(0,0)
	h = np.complex(0,0)
	for j in range(n):
		temp = xh[j]
		h = eps * np.complex(abs(temp.real),abs(temp.imag))
		if abs(h)==0.0:
			h = np.complex(eps, eps)
		xh[j] = temp + h
		h = xh[j] - temp
		f = np.zeros(n,dtype=dataType)
		f = hddbem.bem(xh)
		for i in range(n):
			df[i][j] = (f[i] - fvec[i])/h
		xh[j] = temp
	return df

def fmin(x):
	n=np.shape(x)[0]
	#n=x.shape()[0]
	fvec=hddbem.bem(x)
	if(cplx):
		ssum=complex(0.0,0.0)
	else:
		ssum=0
	for i in range(n):
		ssum+=fvec[i]*fvec[i]
	#print(0.5*np.sqrt(ssum))
	return 0.5*ssum


def lnsrch(xold,fold,g,p,stpmax):
	check=False
	stpmax=100
	ALF = 1.0e-6
	TOLX=1.0e-6
	a=0.0
	alam=0.0
	alam2=0.0
	alamin=0.0
	b=0.0
	dics=0.0
	f2=0.0
	rhs1=0.0
	rhs2=0.0
	slope=0.0
	ssum=0.0
	temp=0.0
	#test=0.0
	tmplam=0.0
	#n = xold.shape()[0]
	n=np.shape(xold)[0]
	for i in range(n):
		ssum+=p[i]*p[i]
	ssum=np.sqrt(ssum)
	if(ssum>stpmax):
		for i in range(n):
			p[i]*=stpmax/ssum
	for i in range(n):
		slope+=g[i]*p[i]
	if(slope>=0.0):
		print("Roundoff error in lnsrch")
	test=0.0
	for i in range(n):
		temp=abs(p[i])/max(abs(xold[i]),1.0)
		if(temp>test):
			test=temp
	alamin=TOLX/test
	alam=1.0

	x=np.zeros(n,dtype=dataType)
	#print("x in lnsrch",x)
	while(True):
		for i in range(n):
			x[i]=xold[i]+alam*p[i]
			f=fmin(x)#fmin
			if(alam<alamin):
				for i in range(n):
					x[i]=xold[i]
				check=True
				#break
				return p,x,f,check
			elif(f<=fold+ALF*alam*slope):
				#break
				return p,x,f,check
			else:
				if(alam==1.0):
					tmplam=-slope/(2.0*(f-fold-slope))
				else:
					rhs1=f-fold-alam*slope
					rhs2=f2-fold-alam2*slope
					a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2)
					b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))
					if(a==0.0):
						tmplam=-slope(2.0*b)
					else:
						disc=b*b-3.0*a*slope
						if(disc<0.0):
							tmplam=0.5*alam2
						elif(b<=0.0):
							tmplam=(-b+np.sqrt(disc))/(3.0*a)
						else:
							tmplam=-slope/(b+np.sqrt(disc))
						if(tmplam>0.5*alam):
							tmplam=0.5*alam
				alam2=alam
				f2=f
				alam=max(tmplam,0.1*alam)	
	##print("x in lnsrch",x)
	return p,x,f,check					



def newt(x):
	check=False
	MAXITS=5000
	TOLF=1.0e-6
	TOLMIN=1.0e-6
	STPMX=100.0
	TOLX=1.0e-6
	n=np.shape(x)[0]
	g=np.zeros(n,dtype=dataType)
	p=np.zeros(n,dtype=dataType)
	xold=np.zeros(n,dtype=dataType)
	fjac=np.zeros((n,n),dtype=dataType)
	#print("#x 1",x)
	if(cplx):
		fjac=fdjac2(x)
	else:
		fjac=fdjac(x)
	fvec=hddbem.bem(x)
	f=fmin(x)
	#print("#f ",f)
	#print("fvec: ",fvec)

	test=0.0
	for i in range(n):
		if(abs(fvec[i])>test):
			test=abs(fvec[i])
	if(test<0.01*TOLF):
		check=False
		return check
	ssum=0.0
	for i in range(n):
		ssum+=x[i]*x[i]
	ssum=np.sqrt(ssum)

	stpmax=STPMX*max(ssum,n)
	#print(ssum)
	#print("stpmax =",stpmax)
	for its in range(MAXITS):
		#print("# its: ", its)
		fvec= hddbem.bem(x)
		if(cplx):
			fjac=fdjac2(x)
		else:
			fjac=fdjac(x)
		for i in range(n):
			ssum=0.0
			for j in range(n):
				ssum+=fjac[j,i]*fvec[j]
			g[i]=ssum
		for i in range(n):
			xold[i]=x[i]
		fold=f
		for i in range(n):
			p[i]=-fvec[i]
		#alu(fjac)
		#alu.solve(p,p)
		#print("p:",p)
		p=np.linalg.solve(fjac,p)
		#print("fjac:",fjac)
		#print("p:",p)
		#print("#x 4",x)
		p,x,f,check=lnsrch(xold,fold,g,p,stpmax)
		#print("x: ",x)
		#print("#x 5",x)
		if(max(abs((x-xold)))<TOLMIN and f<fold):
			break

		print("  ",its+1," ", x[0], "  ",abs((x[0]-xold[0])/xold[0]))
	#if(its+2>MAXITS):
	#	print("MAXITS exceeded in newt")
	#print(x)
	print("===================================================")
	return x
def is_pos_def(mat):
	return np.all(np.linalg.eigvals(mat) > 0)

def newton(x):
	jaco_mat_all=np.zeros((niter,ng+1,ng+1),dtype=dataType)
	x_new = np.ones(ng+1,dtype=dataType)
	for i in range(niter):
		if(cplx):
			jaco_mat = fdjac2(x)
		else:
			jaco_mat = fdjac(x)

		fx = np.transpose(hddbem.bem(x))
		x_new = x - omega * np.dot(np.linalg.inv(jaco_mat),fx)

		k_error = abs((x_new[0].real-x[0].real)/x_new[0].real)
		max_mode_error=max(abs((x_new[1:].real-x[1:].real)/x_new[1:].real))
		max_func_val=max(fx.real)
		
		if(cplx and latex):
			print("  ",str(i+1).zfill(3),"&"," %.5f"%(x_new[0].real),
				"&", "  %.5e"%k_error,"&", "  %.5e"%max_mode_error, 
				"&","  %.5e"%max_func_val,"\\\\")
		elif(cplx and not latex):
			print("  ",str(i+1).zfill(3)," %.5f"%(x_new[0].real), 
				"  %.5e"%k_error, "  %.5e"%max_mode_error, "  %.5e"%max_func_val)
		else:
			print("  ",i+1," ", x_new[0], "  ",abs((x_new[0]-x[0])/x[0]))
		# if ((k_error < eps) and (max_func_val<0.01*eps)):
		if ((k_error < eps) and (max_mode_error<100*eps) and (max_func_val<eps)):
			ctrol.converged=True
			print("="*60)
			break
		x = x_new
		if(i==niter-1):
			print("="*60)
			print("Maximum iteration control is" , niter)
	return x_new,i

def eigen(keff, pplot):
	ctrol.pplot=pplot
	print("===================================================")
	print("         ",ng," group eigenvalue calculation    ")
	#print(" 	                                               ")
	print("="*60)
	print("   iter     keff     err[keff]   max_err[psi]    max_F")
	print("="*60)
	sseed = int(np.loadtxt('.seed')[0])
	if(cplx):
		np.random.seed(sseed)
		x0 = np.ones((nmed-1)*ng+1,dtype=dataType)
		tmp=np.random.random((nmed-1)*ng+1)
		x0.real=tmp/sum(tmp)
		x0.imag=tmp/sum(tmp)
		x0[0] = np.complex(keff,0)
		# print(x0[1:])
#		x0=hddbem.get_initial()
	else:
		np.random.seed(sseed)
		tmp=np.random.random(ng+1)
		x0=tmp/sum(tmp)
		x0[0]=keff
		#print(x0)
	#print(x0)
	#if(ng>25):
	#res=newt(x0)
	#else:
	res,itimes = newton(x0)
	hddbem.save_converged(res)
	if(ctrol.pplot):
		hddbem.get_flx(res)
	return itimes,res[0]
def root_finding(keff):
	print("="*60)
	print("         ",ng," group eigenvalue calculation    ")
	print("="*60)	
	#np.random.seed(4)
	x0 = np.ones(ng+1)
	tmp=np.random.random(ng+1)
	x0=tmp/sum(tmp)
	x0[0] = keff
	# print(x0)

	# method = 
# hybr				not working
# lm                not working
# o broyden1			not working
# o broyden2			sometimes works but not converge
# o anderson			not working
# linearmixing		works but never converge
# diagbroyden 		not working
# excitingmixing 	works but never converge
# o krylov 			not working
# o df-sane 	 	 	works and sometimes converge
	my_options={'line_search':'cruz'}
	# sol = optimize.root(hddbem.bem, x0,method='df-sane',options=my_options)
	sol = optimize.root(hddbem.bem, x0, method='krylov')
	print(sol)



if __name__ == '__main__':
	# inikeff=1.06
	inikeff=1.600879
	if(ng>4):
		iiter,k=eigen(keff=inikeff,pplot=True)
	else:
		#root_finding(keff=inikeff)
		#pass
		iiter,k=eigen(keff=inikeff,pplot=True)
#	print(iiter+1," ",k.real)
	# flx_file="/Users/tian/hddbem/xs/"+str(ng)+"g/flux_bem.txt"
	# cmd='cp .flux.txt '+flx_file
	# os.system(cmd)

