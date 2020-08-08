import numpy as np
import cmath

cplx=True
if(cplx):
	dataType = 'complex'
else:
	dataType = 'float'


def simple_func(z):
	n= np.shape(z)[0]
	I=complex(0,1)
	f=np.zeros(n, dtype=dataType)
	f[0]=np.exp(I*(z[0]+z[1]))
	f[1]=z[0]*z[1]-complex(4,5)
	return f

def function_fx(x):
	size = np.shape(x)[0]
	y = np.zeros(size, dtype=dataType)
	if(cplx):
	    y[0] = x[0] * x[1] * x[2] - complex(6,0);
	    y[1] = x[0] * x[0] * x[1] +  x[1] * x[1] * x[2] +  x[2] * x[2] * x[0] - complex(23,0);
	    y[2] = np.exp(x[0] + x[1] + x[2]) - complex(403,0);
	else:
	    y[0] = x[0] * x[1] * x[2] - 6;
	    y[1] = x[0] * x[0] * x[1] + x[1] * x[1] * x[2] + x[2] * x[2] * x[0] - 23;
	    y[2] = np.exp(x[0] + x[1] + x[2]) - 403;
	return y

def ffdjac(x):
#calculate jacobian matrix for real value function
	EPS = 1.0e-8
	#print(np.shape(x)[0])
	#print("#",x)
	size=len(x)
	#size = np.shape(x)[0]
	xh=np.zeros(size)
	for i in range(size):
		xh[i] = x[i]
	df = np.ones((size,size))
	#print(x)
	fvec = function_fx(x)

	for j in range(size):
		temp = xh[j]
		h = EPS * abs(temp)
		if abs(h)==0.0:
			h = EPS
		xh[j] = temp + h
		#print(xh)
		f = function_fx(xh)
		for i in range(size):
			df[i][j] = (f[i] - fvec[i])/h
			#print(f[i]==function_fx(xh)[i])
			#print(f[i]==function_fx(x)[i])
			#df[i][j] = (function_fx(xh)[i] - function_fx(x)[i])/h
		xh[j] = temp
#	print(f)
#	print(fvec)
	return df


def fdjac2(x):
#calculate jacobian matrix for complex value function
	EPS = 1.0e-8; size = np.shape(x)[0]; xh = x
	df = np.ones((size,size), dtype=complex); 
	#fvec = function_fx(x)
	fvec=simple_func(x)
	for j in range(size):
		temp = xh[j]
		h = EPS * complex(abs(temp.real),abs(temp.imag))
		if abs(h)==0.0:
			h = complex(EPS, EPS)
		xh[j] = temp + h; 	h = xh[j] - temp

		#f = function_fx(xh);
		f=simple_func(xh)
		xh[j] = temp
		for i in range(size):
			df[i][j] = (f[i] - fvec[i])/h
	return df


def fmin(x):
	n=np.shape(x)[0]
	#n=x.shape()[0]
	fvec=function_fx(x)
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
	ALF = 1.0e-9
	TOLX=1.0e-9
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
	if(cplx):
		ssum=complex(0,0)
		temp=complex(0,0)
	else:
		ssum=0
		temp=0
	#test=0.0
	tmplam=0.0
	#n = xold.shape()[0]
	n=np.shape(xold)[0]
	#print("n=",n)
	for i in range(n):
		ssum+=p[i]*p[i]
	ssum=np.sqrt(ssum)
	if(ssum.real>stpmax.real):
		for i in range(n):
			p[i]*=stpmax/ssum
	for i in range(n):
		slope+=g[i]*p[i]
	if(slope>=0.0):
		print("Roundoff error in lnsrch")
	if(cplx):
		test=complex(0,0)
	else:
		test=0.0
	for i in range(n):
		if(cplx):
			temp=complex(abs(p[i].real)/max(abs(xold[i].real),1.0),abs(p[i].imag)/max(abs(xold[i].imag),1.0))
		else:
			temp=abs(p[i])/max(abs(xold[i]),1.0)
		if(abs(temp)>abs(test)):
			test=temp
	if(cplx):
		alamin=TOLX/test.real
	else:
		alamin=TOLX/test
	alam=1.0

	x=np.zeros(n, dtype=dataType)
	#print("x in lnsrch",x)
	while(True):
		for i in range(n):
			#print("x=",x)
			#print("xold=",xold)
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
	MAXITS=2000
	TOLF=1.0e-8
	TOLMIN=1.0e-9
	STPMX=100.0
	TOLX=1.0e-9
	n=np.shape(x)[0]
	g=np.zeros(n, dtype=dataType)
	p=np.zeros(n, dtype=dataType)
	xold=np.zeros(n, dtype=dataType)
	fjac=np.zeros((n,n), dtype=dataType)
	#print("#x 1",x)
	if(cplx):
		fjac=fdjac2(x)
	else:
		fjac=ffdjac(x)
	fvec=function_fx(x)
	#print(x)

	f=fmin(x)
	#print("#f ",f)
	#print("fvec: ",fvec)

	if(cplx):
		test=complex(0.0,0.0)
	else:
		test=0
	for i in range(n):
		if(abs(fvec[i])>test):
			test=abs(fvec[i])
	if(test<0.01*TOLF):
		check=False
		return check

	if(cplx):
		ssum=complex(0.0,0.0)
	else:
		ssum=0
	for i in range(n):
		ssum+=x[i]*x[i]
	ssum=np.sqrt(ssum)
	#ssum=n.sqrt(s)

	stpmax=STPMX*max(ssum.real,n)
	#print(ssum)
	#print("stpmax =",stpmax)
	for its in range(MAXITS):
		#print("# its: ", its)
		fvec= function_fx(x)
		if(cplx):
			fjac=fdjac2(x)
		else:
			fjac=ffdjac(x)
		for i in range(n):
			if(cplx):
				ssum=complex(0.0,0.0)
			else:
				ssum=0
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
		#print(x)
		p,x,f,check=lnsrch(xold,fold,g,p,stpmax)
		print("iter ",its," x: ",x)
		#print("#x 5",x)
		#print(x)
		# if(max(abs((x-xold)))<TOLMIN and f<fold):
		# 	break
		test=0.0
		for i in range(n):
			if(abs(fvec[i])>test):
				test=abs(fvec[i])
		if(test<TOLF):
			check=False
			break
			#return check
		if(check):
			test=0.0
			den=max(f,0.5*n)
			for i in range(n):
				temp=abs(g[i])*max(abs(x[i]),1.0)/den
				if(temp>test):
					test=temp
				check=(test<TOLMIN)
				check=False
				break
				#return check
			test=0.0
			for i in range(n):
				temp=(abs(x[i]-xold[i]))
				if(temp>test):
					test=temp
			if(test<TOLX):
				check=False
				break
				#return check
	if(its+2>MAXITS):
		print("MAXITS exceeded in newt")
	#print(x)
	return x, check	

														
def naive_newt(xx):
	x0=xx
	#x0 = np.transpose(np.array([1,2,-1], dtype=dataType))
	x = x0
	x_new = x0
	# control parameters
	alpha = 1.0
	eps = 1.0e-8
	inter_num = 100
	# root finding with newton's method
	for i in range(inter_num):
		#jacobian_matrix = jacobian(x)
		if(cplx):
			jacobian_matrix=fdjac2(x)
		else:
			jacobian_matrix=ffdjac(x)
		print(jacobian_matrix)
		function_matrix = np.transpose(function_fx(x))
		x_new = x - alpha * np.dot(np.linalg.inv(jacobian_matrix),function_matrix)
		print("iter: ",i," x: ", x_new)
		error = max(abs(x_new-x))
		#print("error: ", error)
		if error < eps:
			print("Totoal interation times: ", i+1); break
		x = x_new
		
	return x
			#print(function_fx(x_new))

def test_newton():
	size=3
	xx=np.ones(size, dtype=dataType)
	#np.random.seed(3)
	for i in range(size):
		if(cplx):
			xx[i]=complex(i+1)
		else:
			xx[i]=i+1
	x2=naive_newt(xx)
	print("Solution by Naive newton: ",x2)
	print(function_fx(x2))

	# x1=newt(xx)
	# print("Solution by Newton: ",x1[0])
	# print(function_fx(x1[0]))
def test_simplex_complex_funtion():
	z=np.zeros(2, dtype=dataType)
	z[0]=complex(1,1)
	z[1]=complex(-1,-2)
	
	jaco=fdjac2(z)
	print(jaco)

test_simplex_complex_funtion()
