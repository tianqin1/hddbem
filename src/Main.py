import numpy as np
import cmath
import Control as ctrol
from Control import *
import HDDBEM as hddbem
#from TwoGroup import *


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
# this is line 40
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
		# test=0.0
		# for i in range(n):
		# 	if(abs(fvec[i])>test):
		# 		test=abs(fvec[i])
		# if(test<TOLF):
		# 	check=False
		# 	break
		# 	#return check
		# if(check):
		# 	test=0.0
		# 	den=max(f,0.5*n)
		# 	for i in range(n):
		# 		temp=abs(g[i])*max(abs(x[i]),1.0)/den
		# 		if(temp>test):
		# 			test=temp
		# 		check=(test<TOLMIN)
		# 		check=False
		# 		break
		# 		#return check
		# 	test=0.0
		# 	for i in range(n):
		# 		temp=(abs(x[i]-xold[i]))
		# 		if(temp>test):
		# 			test=temp
		# 	if(test<TOLX):
		# 		check=False
		# 		break
		# 		#return check
		print("  ",its+1," ", x[0], "  ",abs((x[0]-xold[0])/xold[0]))
	#if(its+2>MAXITS):
	#	print("MAXITS exceeded in newt")
	#print(x)
	print("===================================================")
	return x


def newton(x):
	x_new = np.ones(ng+1,dtype=dataType)
	for i in range(niter):
		if(cplx):
			jaco_mat = fdjac2(x)
		else:
			jaco_mat = fdjac(x)
		#print(jaco_mat)
		fx = np.transpose(hddbem.bem(x))
		x_new=x-np.linalg.solve(jaco_mat,fx)
		#x_new = x - alpha * np.dot(np.linalg.inv(jaco_mat),fx)
		# if(cplx):
		# 	print("  ",i+1," ", x_new[0].real, "  ",
		# abs((x_new[0].real-x[0].real)/x[0].real))
		# else:
		# 	print("  ",i+1," ", x_new[0], "  ",abs((x_new[0]-x[0])/x[0]))
		if(cplx):
			print("  ",i+1," ", x_new.real, "  ",
					abs((x_new.real-x.real)/x.real))
		else:
			print("  ",i+1," ", x_new[0], "  ",
				abs((x_new[0]-x[0])/x[0]))

		error = max(abs((x_new-x)/x))
		#print("error: ", error)
		if error < eps:
			ctrol.converged = True
			#print(converged)
			print("===================================================")
			#print("Totoal interation times: ", i+1)
			#print("===================================================")
			break
		x = x_new
	return x_new

def eigen(keff, pplot):
	ctrol.pplot=pplot
	#print("===================================================")
	print("         ",ng," group eigenvalue calcaution    ")
	#print(" 	                                               ")
	print("===================================================")
	print("  iter        keff                 err[keff]       ")
	print("===================================================")
	if(cplx and ng<=10):
		#print("==>> Run at complex mode")
		x0 = np.ones(ng+1,dtype=dataType)
		x0[0] = np.complex(keff,0)
	elif(cplx and ng>10):
		x0 = 1/ng*np.ones(ng+1,dtype=dataType)
		x0.real=1/ng*np.random.random(ng+1)
		x0.imag=1/ng*np.random.random(ng+1)
		x0[0] = np.complex(keff,0)
	elif(ng==107):
		inmode=hddbem.get_inmode(ng,keff)
		x0=np.ones(ng+1,dtype=dataType)
		for i in range(1,ng+1):
			x0[i]=inmode[i-1]
		x0[0] = np.complex(keff,0)
	elif(not cplx and ng<10):
		x0 = 0.05*np.ones(ng+1)
		x0[0] = keff
	else:
		x0=np.random.random(ng+1)
		x0[0]=keff
	#print(x0)
	#if(ng>25):
	#res=newt(x0)
	#else:
	res = newton(x0)
	if(ctrol.pplot):
		hddbem.get_flx(res)

if __name__ == '__main__':
	inikeff=1.0
	if(ng>2 and ng < 30):
		inikeff=1.22268727
	elif(ng==107):
		inikeff=1.21763 
	else:
		inikeff=1.35
	#eigen(keff=inikeff,pplot=False)
	eigen(keff=inikeff,pplot=True)

