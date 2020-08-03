import math
from pylab import *
import scipy.integrate as integrate
def epower(p):
	y=[]
	for i in range(len(p)):
		y.append(exp(p[i]))               # exp(x) function
	return y
def coscos(p):
	z=[]
	for i in range(len(p)):               # cos(cos(x))
		z.append(cos(cos(p[i])))
	return z
t=linspace(-2*(pi),4*(pi),401)
t=t[:-1]
ft=linspace(0,2*(pi),401)
ft=ft[:-1]
no_of_coefficients=25
c=1/math.pi
u1=lambda x,k: (exp(x))*cos(k*x)                    # creation of a function u(x,k) to find an coefficients for exp(x)
v1=lambda x,k: (exp(x))*sin(k*x)					# creation of a function v(x,k) to find bn coefficients for exp(x)
u2=lambda x,k: (cos(cos(x)))*cos(k*x)				# creation of a function u(x,k) to find an coefficients for cos(cos(x))
v2=lambda x,k: (cos(cos(x)))*sin(k*x)				# creation of a function v(x,k) to find bn coefficients for cos(cos(x))
a1=[(c/2)*integrate.quad(u1,0,2/c,args=(0,))[0]]
b1=[(c/2)*integrate.quad(v1,0,2/c,args=(0,))[0]]
a2=[(c/2)*integrate.quad(u2,0,2/c,args=(0,))[0]]
b2=[(c/2)*integrate.quad(v2,0,2/c,args=(0,))[0]]
for i in range(1,no_of_coefficients+1):
    a1.append(c*integrate.quad(u1,0,2/c,args=(i,))[0])      #storing an coefficients of exp(x) in a1
    a2.append(c*integrate.quad(u2,0,2/c,args=(i,))[0])		#storing an coefficients of cos(cos(x)) in a2	
    b1.append(c*integrate.quad(v1,0,2/c,args=(i,))[0])		#storing bn coefficients of exp(x) in a1
    b2.append(c*integrate.quad(v2,0,2/c,args=(i,))[0])		#storing bn coefficients of cos(cos(x)) in a2
ans1=[a1[0]]
ans2=[a2[0]]
for i in range(1,len(a1)):
    ans1.append(a1[i])
    ans1.append(b1[i])									#gettting a column matrix of all coefficients of exp(x)
    ans2.append(a2[i])
    ans2.append(b2[i])									#getting a column matrix of all coefficients of cos(cos(x))
p=linspace(0,2*pi,401) 
p=p[:-1]
B1=epower(p)
B2=coscos(p)
A=zeros((400,2*no_of_coefficients+1))
A[:,0]=1
for k in range(1,no_of_coefficients+1):
    A[:,2*k-1]=cos(k*p)
    A[:,2*k]=sin(k*p)									# getting matrix A to find bestfit fourier coefficients
c1=lstsq(A,c_[B1],rcond=1)[0]
c2=lstsq(A,c_[B2],rcond=1)[0]							#getting bestfit coefficients
besta1=[c1[0]]
bestb1=[c1[0]]
besta2=[c2[0]]
bestb2=[c2[0]]
for k in range(1,no_of_coefficients+1):
    besta1.append(c1[2*k-1])							#storing bestfit an coefficients of exp(x) in besta1 
    besta2.append(c2[2*k-1])							#storing bestfit an coefficients of cos(cos(x)) in besta2
    bestb1.append(c1[2*k])								#storing bestfit an coefficients of exp(x) in bestb1
    bestb2.append(c2[2*k])								#storing bestfit an coefficients of cos(cos(x)) in bestb2
gsize=3
figure(3)
semilogy(abs(array(a1)),'ro',label='integrated coeff an')
semilogy(abs(array(b1)),'r+',label='integrated coeff bn')                 #semilogy plot of an and bn coefficients of exp(x)
semilogy(abs(array(besta1)),'go',markersize=gsize,label='bestfit an')	  
semilogy(abs(array(bestb1)),'g+',markersize=gsize,label='bestfit bn')		#semilogy plot of bestfit an and bn coefficients of exp(x)
legend(loc='upper right')
figure(4)
loglog(abs(array(a1)),'ro',label='integrated coeff an')
loglog(abs(array(b1)),'r+',label='integrated coeff bn')						 #loglog plot of an and bn coefficients of exp(x)
loglog(abs(array(besta1)),'go',markersize=gsize,label='bestfit an')
loglog(abs(array(bestb1)),'g+',markersize=gsize,label='bestfit bn')				#loglog plot of bestfit an and bn coefficients of exp(x)
legend(loc='upper right')
figure(5)
semilogy(abs(array(a2)),'ro',label='integrated coeff an')
semilogy(abs(array(b2)),'r+',label='integrated coeff bn')					#semilogy plot of an and bn coefficients of cos(cos(x))
semilogy(abs(array(besta2)),'go',markersize=gsize,label='bestfit an')
semilogy(abs(array(bestb2)),'g+',markersize=gsize,label='bestfit bn')		#semilogy plot of bestfit an and bn coefficients of cos(cos(x))
legend(loc='upper right')
figure(6)
loglog(abs(array(a2)),'ro',label='integrated coeff an')
loglog(abs(array(b2)),'r+',label='integrated coeff bn')						#loglog plot of an and bn coefficients of cos(cos(x))
loglog(abs(array(besta2)),'go',markersize=gsize,label='bestfit an')
loglog(abs(array(bestb2)),'g+',markersize=gsize,label='bestfit bn')			#loglog plot of bestfit an and bn coefficients of cos(cos(x))
legend(loc='upper right')
	
absdiff1=abs(c_[ans1]-c1)
absdiff2=abs(c_[ans2]-c2)
print('for exp(x) largest deviation '+str(max(absdiff1)))					# value of max deviation in fourier coefficients for exp(x)
print('for cos(cos(x)) largest deviation '+str(max(absdiff2)))				# value of max deviation in fourier coefficients for cos(cos(x))
figure(1)
semilogy(t,e**t,color='blue',label='exp(x)')
semilogy(t,epower(abs(t%(2*pi))),color='red',label='expected from fourier')
semilogy(p,dot(A,c1),'go',label='exp(x) from bestfit coeff')				# plotting of function values of exp(x) with bestfit fourier coefficients
xlabel('x ',size=20)
ylabel(r'$exp(x)$',size=20)
legend(loc='upper right')
grid(True)
figure(2)
plot(t,coscos(t),label='cos(cos(x))')
plot(p,dot(A,c2),'go',label='cos(cos(x)) from bestfit coeff')				# plotting of function values of cos(cos(x)) with bestfit fourier coefficients
xlabel('x ',size=20)
ylabel(r'$cos(cos(x))$',size=20)
legend(loc='upper right')
grid(True)
show()