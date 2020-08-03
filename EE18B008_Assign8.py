from pylab import * 

tm=linspace(-4*pi,4*pi,513)
tm=tm[:-1]
w=linspace(-64,64,513)
w=w[:-1]

y=[sin(5*tm),(1+0.1*cos(tm))*cos(10*tm),sin(tm)**3,cos(tm)**3]
yn=['sin(5t)','(1+0.1cos(t))cos(10t)','sin(t)^3','cos(t)^3']
for i in range(len(y)):
	Y=fftshift(fft(y[i]))/512.0 
	figure(i)
	subplot(2,1,1) 
	plot(w,abs(Y),lw=2) 
	xlim([-15,15]) 
	ylabel(r"$|Y|$",size=16) 
	title(r"Spectrum of "+ yn[i]) 
	grid(True) 
	subplot(2,1,2) 
	plot(w,angle(Y),'bo',markersize=2,lw=2) 
	ii=where(abs(Y)>1e-3) 
	plot(w[ii],angle(Y[ii]),'ro',lw=2) 
	xlim([-15,15]) 
	ylabel(r"Phase of $Y$",size=16) 
	xlabel(r"$\omega$",size=16) 
	grid(True)
    
y3=cos(20*tm+5*cos(tm))
Y=fftshift(fft(y3))/512.0
figure(len(y)) 
subplot(2,1,1) 
plot(w,abs(Y),lw=2) 
xlim([-32,32]) 
ylabel(r"$|Y|$",size=16) 
title(r"Spectrum of "+ 'cos(20t+5cos(t))') 
grid(True) 
subplot(2,1,2) 
#plot(w,angle(Y),'ro',lw=2) 
ii=where(abs(Y)>1e-3) 
plot(w[ii],angle(Y[ii]),'ro',lw=2)
xlim([-32,32])  
ylabel(r"Phase of $Y$",size=16) 
xlabel(r"$\omega$",size=16) 
grid(True)

y4=e**(-tm*tm/2)
Y=fftshift(fft(y4))*(8*pi)/(512)    # Multiply by T/N i.e., time period/no. of samples
figure() 
subplot(2,1,1)
plot(w,abs(Y),lw=2) 
xlim([-15,15]) 
ylabel(r"$|Y|$",size=16) 
title(r"Spectrum of "+'e^(-t^2/2)') 
grid(True)  
subplot(2,1,2) 
#plot(w,angle(Y),'ro',lw=2) 
ii=where(abs(Y)>1e-3) 
plot(w[ii],angle(Y[ii]),'ro',lw=2)
xlim([-15,15]) 
ylabel(r"Phase of $Y$",size=16) 
xlabel(r"$\omega$",size=16) 
grid(True)

show()