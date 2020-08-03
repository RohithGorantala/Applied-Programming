import scipy.signal as sp
from pylab import *
c2=0.05                         #decay constant is 0.05
c1=0.5							#decay constant is 0.5
p1=poly1d([1,c1])
p2=poly1d([1,2*c1,2.25+c1**2])
p3=poly1d([1,0,2.25])
d1=poly1d([1,c2])
d2=poly1d([1,2*c2,2.25+c2**2])
d3=poly1d([1,0,2.25])
X1=sp.lti(p1,polymul(p2,p3))        # Transfer function with f(t) having decay constant 0.5
X2=sp.lti(d1,polymul(d2,d3))		# Transfer function with f(t) having decay constant 0.05
t1,x1=sp.impulse(X1,None,linspace(0,50,2000))               	# sp.impulse gives time response x(t)
t2,x2=sp.impulse(X2,None,linspace(0,50,2000))
figure(1)
title('Response of a spring satisfying x" +2.25x=f(t)')
plot(t1,x1,label='at decay const=0.5')
plot(t2,x2,label='at decay const=0.05')
legend(loc='upper left')
w=linspace(1.4,1.6,5)									# taking omega from 1.4 to 1.6 in steps of 0.05
t=t2
y=[]
figure(2)
title('Response of a spring satisfying x" +2.25x=f(t) at different frequencies')
for i in w:
    H=sp.lti([1],[1,0,2.25])
    f=cos(i*t)*(e**(-0.05*t))
    y.append(sp.lsim(H,f,t)[1])							# sp.lsim will simulate the output signal for given transfer function and input
    plot(t,y[-1],label='at w='+str(i))
legend(loc='upper left')

Xs=sp.lti([1,0,2],[1,0,3,0])					# laplace transform of x(t)
Ys=sp.lti([2],[1,0,3,0])						# laplace transform of y(t)
tx,x=sp.impulse(Xs,None,linspace(0,20,2000))
ty,y=sp.impulse(Ys,None,linspace(0,20,2000))
figure(3)
title('Responses x(t) and y(t))')
plot(tx,x,label='x(t)')
plot(ty,y,label='y(t)')
legend(loc='upper left')
'''for the given circuit with values R=100,L=1e-6,c=1e-6, We calculate the transfer function of the circuit'''
Hs=sp.lti([1],[1e-12,1e-4,1]) 
freq,S,phi=Hs.bode()						# Hs.bode will return frequency vector,magnitude vector and phase vector.
figure(4)
subplot(1,2,1)
title('magnitude plot of H(s)')
semilogx(freq,S)							# Bode plot of magnitude of transfer function
subplot(1,2,2)
title('phase plot of H(s)')
semilogx(freq,phi)							# Bode plot of phase of transfer function
figure(5)
'''Listing time vectors with diffrent time steps'''
timee=[linspace(0,3*(10**-5),1000),linspace(0,3*(10**-5),3),linspace(0,10*(10**-3),60),linspace(0,10*(10**-3),120000)]
vo=[]
for i in timee:
    vi=cos((10**3)*i)-cos((10**6)*i)						# input to the circuit 
    vo.append(sp.lsim(Hs,vi,i)[1])							# sp.lsim will simulate the output signal for given transfer function and input
    plot(i,vo[-1],label='time step='+str(i[-1]/len(i)))		# Plotting the output signal
legend(loc='upper left')
figure(7)
tm=linspace(0,12*(10**-3),50000)
vi=cos((10**3)*tm)-cos((10**6)*tm)
plot(tm,sp.lsim(Hs,vi,tm)[1])
figure(6)												# Observing the behaviour of the output signal with different time steps and time scales
subplot(7,1,1)
plot(timee[0],vo[0],label='usec: time step< 1e-6')
legend(loc='upper left')
subplot(7,1,3)
plot(timee[1],vo[1],label='usec: time step> 1e-6')
legend(loc='upper left')
subplot(7,1,5)
plot(timee[2],vo[2],label='msec: time step> 1e-6')
legend(loc='upper left')
subplot(7,1,7)
plot(timee[3],vo[3],label='msec: time step< 1e-6')
legend(loc='upper left')

show()