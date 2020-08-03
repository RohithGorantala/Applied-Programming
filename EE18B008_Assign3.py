from pylab import *
import scipy.special as sp
from scipy.linalg import lstsq 

sigma=logspace(-1,-3,9)
li=loadtxt('fitting.dat')
f1=[[0 for i in range(len(li))] for i in range(10)]
for i in range(10):														# Loading text from fitting.dat and storing each column in list
	for j in range(len(li)):											# Where first column represents time values and other 9 columns 
		f1[i][j]=li[j][i]												# are functions with different noise levels
t=linspace(0,10,len(li))												#time array

figure(1)
for i in range(1,10):
    plot(f1[0],f1[i],label='c'+str(i)+'='+str(round(sigma[i-1],3)))		# Plotting each function with different noise in a plot
    legend(loc='upper right')
a=1.05
b=-0.105
y=a*sp.jn(2,t)+b*t 														# Actual function with no noise
plot(t,y,color='black',linewidth=2,label='true value')
legend(loc='upper right')												# getting legend
xlabel(r'$t -->$',size=20) 
ylabel(r'$f(t)+noise$',size=20)
title(r'Data to be fitted to theory')
grid(True)

figure(2)
i= int(input('give the stdev'))
errorbar(t[::5],f1[1][::5],sigma[i],fmt='ro',label='error bar')			# Plotting Errorbar 
plot(t,y,color='black',linewidth=2,label='f(t)')
legend(loc='upper right')
xlabel(r'$t$',size=20)
title(r'Data points for c='+str(sigma[i])+' along with exact function')
grid(True)

M=array([[sp.jn(2,i),i] for i in t])
p=array(c_[[a,b]])
d=dot(M,p)																# Verifying multiplication of vectors with actual funtion vector
y=array([[y[i]] for i in range(len(y))])
if(d.all()==y.all()):
	print('They are EQUAL')

figure(3)
A=linspace(0,2,21)										# A array
B=linspace(-0.2,0,21)									# B array
e=[[] for i in range(len(A))]
										
for i in range(len(A)):
    for j in range(len(B)):
        v=0
        for k in range(len(t)):
            v=v+pow((f1[1][k]-(A[i]*sp.jn(2,t[k])+B[j]*t[k])),2)
        v=v/101
        e[i].append(v)

A,B=meshgrid(A,B)
c=contour(A,B,array(e),20)								# Contour Plot
clabel(c,c.levels[:5],inline=True,fontsize=10,)				# Labeling of contour levels
plot(a,b,color='red',marker='o')												
text(a,b,'Exact location')
xlabel(r'$A  -->$',size=20)
ylabel(r'$B  -->$',size=20)
title(r'contour plot of Eij')

figure(4)
p=[]
for i in range(9):
    p.append(lstsq(M,f1[i+1]))								# getting best estimate of A and B values
Aerr=[]
Berr=[]
for i in range(len(p)):
    Aerr.append(abs(p[i][0][0]-1.05))
    Berr.append(abs(p[i][0][1]+0.105))
plot(sigma,Aerr,'r--o',linewidth='0.8',label='Aerr')		# Plotting Aerr vs corresponing noise
plot(sigma,Berr,'g--o',linewidth='1',label='Berr')			# Plotting Berr vs corresponiing noise
legend(loc='upper left')
xlabel(r'$Noise Standard deviation -->$',size=10)
ylabel(r'$MS error  -->$',size=10)
title(r'Variation of error with Noise')
grid(True)

figure(5)
loglog(sigma,Aerr,'o')										# Making x and y axes in to logX and logY axes
loglog(sigma,Berr,'o')
errorbar(sigma,Aerr,sigma[0],fmt='ro',label='Aerr')			# Plotting Errorbar
errorbar(sigma,Berr,sigma[0],fmt='go',label='Berr')
legend(loc='upper right')
xlabel(r'$Cn -->$',size=10)
ylabel(r'$MS error  -->$',size=10)
title(r'Variation of error with noise')
grid(True)
show()