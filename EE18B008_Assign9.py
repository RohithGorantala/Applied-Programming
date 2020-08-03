from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3

def get_data(w,Y,titlee,xlimit):					# Function to plot DFT of the signal
    figure()
    subplot(211)
    plot(w,abs(Y),lw=2)
    title(titlee)
    ylabel(r"$|Y|$",size=16)
    xlim([-xlimit,xlimit])
    grid(True)
    subplot(212)
    ii=where(abs(Y)>1e-3)
    '''  if(plotgo==1):
        plot(w[ii],angle(Y[ii]),'go',lw=2)
    else:
        plot(w,angle(Y),'ro',lw=2) '''
    plot(w,angle(Y),'ro',lw=2)
    xlabel(r"$omega$",size=16)
    ylabel(r"Phase of $Y$",size=16)
    xlim([-xlimit,xlimit])
    grid(True)
    show()

def fft_wnd(t,y,window_flag):							# Function which take cares of windowing
    n=arange(len(t))
    wnd=fftshift(0.54+0.46*cos(2*pi*n/(len(t)-1)))
    if(window_flag==1):
        y=y*wnd
        y[0]=0
        y=fftshift(y)
        Y=fftshift(fft(y))/len(y)
        return Y
    else:
        y[0]=0
        y=fftshift(y)
        Y=fftshift(fft(y))/len(y)
        return Y

def estimator(t,w,omega,delta,window_flag,noise_flag):	# Function which estimates where the peaks are and also error between actual and calculated values.
    if(noise_flag==1):
        y=cos(omega*t+delta)+0.1*randn(len(t3))
    else:
        y=cos(omega*t+delta)
    Y=fft_wnd(t,y,window_flag)
    ii = np.where(w>0)[0]
    ii = np.where((abs(Y) == max(abs(Y[ii]))))
    est_delta = abs(np.angle(Y[ii])[0])   # The phase of the graph at the peak is the delta
    ii = np.where((abs(Y) > 3.5e-2) & (w >= 0))[0]
    est_omega = abs(Y[ii]*w[ii])
    est_omega = sum(est_omega)/(sum(abs(Y[ii])))# As peak is spread out, omega is estimated as the weighted average(centre of mass) of the broad area near the peak.
    get_data(w,Y,'Spectrum of $cos($'+str(omega)+'$t$'+'+'+str(delta)+'$)$',omega+4)
    print ('the Calculated delta is %.6f and the error in the calculated delta is %.6f' %(est_delta, abs(est_delta - delta)))
    print ('the Calculated omega is %.6f and the error in the calculated omega is %.6f' %(est_omega, abs(est_omega - omega)))


# Problem 1 and 2

t1=linspace(-pi,pi,65)[:-1]
fmax1=1/(t1[1]-t1[0])
t2=linspace(-4*pi,4*pi,257)[:-1]
fmax2=1/(t2[1]-t2[0])
w1=linspace(-pi*fmax1,pi*fmax1,65)[:-1]
w2=linspace(-pi*fmax2,pi*fmax2,257)[:-1]

Y=fft_wnd(t2,cos(0.86*t2)**3,0)									# DFT of cos(0.86t)^3 without window
get_data(w2,Y,r"Spectrum of cos(0.86t)^3 without window",5)		# plots DFT of given Y

Y=fft_wnd(t2,cos(0.86*t2)**3,1)									# DFT of cos(0.86t)^3 with window
get_data(w2,Y,r"Spectrum of cos(0.86t)^3 with window",5)		# plots DFT of given Y

# Problem 3 and 4

t3=linspace(-pi,pi,129)[:-1]
fmax3=1/(t3[1]-t3[0])
w3=linspace(-pi*fmax3,pi*fmax3,len(t3)+1)[:-1]

omega=float(input("Enter your value of omega between 0.5 and 1.5"))				# Takes input for omega
delta=float(input("Enter your value of delta"))									# Takes input for delta
noise_flag=float(input("Enter 1 to add noise or else 0"))						# Takes input for noise addition

estimator(t3,w3,omega,delta,1,noise_flag)# plots DFT for given omega and delta and estimates where the peaks are and also error in estimated and given values.

# Problem 5
t4=linspace(-pi,pi,1025)[:-1]
fmax4=1/(t4[1]-t4[0])
w4=linspace(-pi*fmax4,pi*fmax4,len(t4)+1)[:-1]
y=cos(16*(1.5+(t4/(2*pi)))*t4)

Y=fft_wnd(t4,y,1)												# gets DFT of chirped signal
get_data(w4,Y,'spectrum of $cos(16(1.5+t/2*pi))t$',100)			# plots DFT of chirped signal

# Problem 6

w5=linspace(-pi*fmax4,pi*fmax4,65)[:-1]
t5=linspace(-pi,pi,17)[:-1]
tt,ww=meshgrid(t5,w5)
ns=arange(64)
wnd=fftshift(0.54+0.46*cos(2*pi*ns/63))

X=[fftshift(fft(fftshift(y[0:64]*wnd)))/64]
for i in range(1,16):
    X.append(fftshift(fft(fftshift(y[i*64:(i+1)*64]*wnd)))/64)			# finding DFT's for each 64 samples of signal

Y=X[0]
for i in range(1,16):
    Y=c_[Y,X[i]]														# Storing those 16 DFT sets as columns in 2D-array

figure(6)
ax=p3.Axes3D(figure(6))
title("3D-surface plot of DFT with window")
ax.set_xlim3d(-80,80)
ax.set_xlabel(r"frequency")
ax.set_ylabel(r"time")
ax.set_zlabel(r"DFT Spectrum")
ww[ww>80]= nan      # Without this and the next line, the surface plot overflows due to the setting of xlim.
ww[ww<-80]= nan
surf = ax.plot_surface(ww,tt,abs(Y),rstride=1,cstride=1,cmap=cm.jet,linewidth=0,antialiased=False)# 3D surface plot of splitted DFTs vs time and frequency.
show()