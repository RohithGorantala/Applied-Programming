from sympy import *
import scipy.signal as sp
import pylab as p

s=symbols('s')

def lowpass(R1,R2,C1,C2,G,Vi):	                     # Function to get transfer function of lowpass circuit

	A=Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
	b=Matrix([0,0,0,Vi/R1]) 
	V=A.inv()*b
	return (A,b,V)
A,b,V=lowpass(10000,10000,1e-9,1e-9,1.586,1)
Vo=V[3]
ww=p.logspace(0,8,801) 
ss=1j*ww
hf=lambdify(s,Vo,'numpy')                           # makes it as python function
v=hf(ss)
h2=simplify(Vo)                                     # Simplifies it and returns reduced form
numer,denom=h2.as_numer_denom()						# Extracts numerator and denominator
num=poly(numer,s)
den=poly(denom,s)									# making them as polynomials to extract coefficients using all_coeff()

Hlp=sp.lti([float(i) for i in num.all_coeffs()],[float(i) for i in den.all_coeffs()])      # Impulse response of Lowpass filter
p.figure(1)
p.title("impulse reponse of lowpass filter")
p.xlabel(r'$s$',size=10)
p.loglog(ww,abs(v),lw=2)							# Frequency response of lowpass filter
p.grid(True)

def highpass(R1,R3,C1,C2,G,Vi):                      # function to get transfer function of highpass circuit

    A=Matrix([[0,0,1,-1/G],[0,G,-G,-1],[s*C2*R3/(1+s*C2*R3),-1,0,0],[(1/R1)+s*C1,1/R3,0,-1/R1]])
    b=Matrix([0,0,0,Vi*s*C1])
    V=A.inv()*b
    return (A,b,V)
A,b,V=highpass(10000,10000,1e-9,1e-9,1.586,1)
Vo=V[3]
hf=lambdify(s,Vo,'numpy') 
v=hf(ss)
h2=simplify(Vo)
numer,denom=h2.as_numer_denom()
num=poly(numer,s)
den=poly(denom,s)
Hhp=sp.lti([float(i) for i in num.all_coeffs()],[float(i) for i in den.all_coeffs()])      # Impulse response of Highpass filter
p.figure(2)
p.title("impulse response of Highpass filter")
p.xlabel(r'$s$',size=10)
p.loglog(ww,abs(v),lw=2)							# Frequency response of highpass filter
p.grid(True)

tm=p.linspace(0,1.5e-3,50000)

v1=p.array([1 for i in tm])						# Step function
y1=sp.lsim(Hlp,v1,tm)[1]						# simulate output for given input and system transfer function
print(y1)
p.figure(3)
p.title("step reponse of lowpass filter")
p.xlabel(r'$t$',size=10)
p.plot(tm,y1)

v2=p.sin(2000*p.pi*tm)+p.cos(2*1e6*p.pi*tm)     # given input
y2=sp.lsim(Hlp,v2,tm)[1]
p.figure(4)
p.title("input reponse of lowpass filter")
p.xlabel(r'$t$',size=10)
p.plot(tm,y2)

y3=sp.lsim(Hhp,v2,tm)[1]
p.figure(5)
p.title("input response of highpass filter")
p.xlabel(r'$t$',size=10)
p.plot(tm,y3)

omega=500
decay=0.05
v4=(p.e**(-decay*tm))*p.sin(omega*tm)			# Damped sinusoid
y4=sp.lsim(Hhp,v4,tm)[1]
p.figure(6)
p.title("response of highpass filter for damped sinusoid")
p.xlabel(r'$t$',size=10)
p.plot(tm,y4)

y5=sp.lsim(Hhp,v1,tm)[1]					
p.figure(7)
p.title("step response of highpass filter")
p.xlabel(r'$t$',size=10)
p.plot(tm,y5)

p.show()