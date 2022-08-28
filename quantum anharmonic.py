import numpy as np
import matplotlib.pyplot as p
Ein=.1
Ef=22.0
L=10
k=.5
def W(e,x):
	return e-(0.5*k*x**2+b*x**4+c*x**6)
def phiEo(E):
	psi1=0
	psi2=1.0
	sum=0
	h=0.001
	xc=np.sqrt((2*abs(E))/k)
	n=int((xc+L)/h+0.5)
	for i in range(0,n+1):
		psi3=(2*psi2-psi1-(h**2/12.0)*(10*W(E,-L+(1+i)*h)*psi2+W(E,-L+i*h)*psi1))/(1.0+(h**2/12.0)*W(E,-L+(2+i)*h))
		psi1=psi2
		psi2=psi3
	#print(psi1)
	o1=psi1
	p1=psi2
	psi1=0
	psi2=1.0
	m=int((L-xc)/h+0.5)
	for i in range(0,m+1):
		psi3=(2*psi2-psi1-(h**2/12.0)*(10*W(E,L-(1+i)*h)*psi2+W(E,L-i*h)*psi1))/(1.0+(h**2/12.0)*W(E,L-(2+i)*h))
		psi1=psi2
		psi2=psi3
	#print(psi1)
	o2=psi1
	p2=psi2
	#print(p1,p2)
	s=o1/o2
	#print (p2-s*p2)
	return p2-s*p1
def sec(f,ig):
	x1=ig
	x2=ig+0.1
	f1 = f(x1)
	f2 = f(x2)
	#print(f1,f2)
	#print(x1,x2)
	x3= (f2*x1 - f1*x2)/(f2 - f1)
	#print(x3)
	while abs(x3-x2)>=0.0001:
		x1=x2
		x2=x3
		f1 = f(x1)
		f2 = f(x2)
		x3= (f2*x1 - f1*x2)/(f2 - f1)
		#print(f1,f2)
		#print(x1,x2)
		#print(x3)
	return x3
def norm(E):
	sum=0
	psi1=0
	psi2=1.0
	sum=0
	h=0.001
	xc=np.sqrt((2*abs(E))/k)
	n=int((xc+L)/h+0.5)
	for i in range(0,n+1):
		sum=sum+psi1*psi1
		psi3=(2*psi2-psi1-(h**2/12.0)*(10*W(E,-L+(1+i)*h)*psi2+W(E,-L+i*h)*psi1))/(1.0+(h**2/12.0)*W(E,-L+(2+i)*h))
		psi1=psi2
		psi2=psi3
	#print(psi1)
	o1=psi1
	p1=psi2
	psi1=0
	psi2=1.0
	m=int((L-xc)/h+0.5)
	for i in range(0,m+1):
		sum=sum+psi1*psi1
		psi3=(2*psi2-psi1-(h**2/12.0)*(10*W(E,L-(1+i)*h)*psi2+W(E,L-i*h)*psi1))/(1.0+(h**2/12.0)*W(E,L-(2+i)*h))
		psi1=psi2
		psi2=psi3
	#print(psi1)
	o2=psi1
	p2=psi2
	#print(p1,p2)
	s=o1/o2
	#print (p2-s*p2)
	sum=1/np.sqrt(sum*h)
	return s,sum
def phi(E):
	s,n1=norm(E)
	sum=0
	psi1=0
	psi2=1.0
	sum=0
	h=0.001
	psi=[]
	x=[]
	xc=np.sqrt((2*E)/k)
	n=int((xc+L)/h+0.5)
	for i in range(0,n+1):
		sum=sum+psi1*psi1
		psi3=(2*psi2-psi1-(h**2/12.0)*(10*W(E,-L+(1+i)*h)*psi2+W(E,-L+i*h)*psi1))/(1.0+(h**2/12.0)*W(E,-L+(2+i)*h))
		psi1=psi2
		psi2=psi3
		psi.append(psi1*n1)
		x.append(-L+i*h)
	#print(psi1)
	o1=psi1
	p1=psi2
	psi1=0
	psi2=1.0
	m=int((L-xc)/h+0.5)
	for i in range(0,m+1):
		sum=sum+psi1*psi1
		psi3=(2*psi2-psi1-(h**2/12.0)*(10*W(E,L-(1+i)*h)*psi2+W(E,L-i*h)*psi1))/(1.0+(h**2/12.0)*W(E,L-(2+i)*h))
		psi1=psi2
		psi2=psi3
		psi.insert(n+1,psi1*n1*s)
		x.insert(n+1,L-i*h)
	#print(psi1)
	o2=psi1
	p2=psi2
	#print(p1,p2)
	s=o1/o2
	#print (p2-s*p2)
	p.plot(x,psi)
b=0.0
c=0.0
En=sec(phiEo,.5)
print(En)
phi(En)
b=0.2
c=0.0
En=sec(phiEo,.5)
print(En)
phi(En)
b=.2
c=0.009
En=sec(phiEo,.5)
print(En)
phi(En)
p.xlabel("x")
p.ylabel("psi")
p.legend(("Harmonic","Quartic","sextic"))
p.title("Wave function of anharmonic QHO")
p.savefig("QHO anharmonic")
p.show()

