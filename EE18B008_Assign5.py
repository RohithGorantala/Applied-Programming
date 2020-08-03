from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import os,sys

if(len(sys.argv)==5):
	Nx=int(sys.argv[1])
	Ny=int(sys.argv[2])
	if(int(sys.argv[3])<((int(sys.argv[2])-1)/2)+1):
		radius=int(sys.argv[3])
	else:
		print("wrong input")
		exit()
	Niter=int(sys.argv[4])
else:
	print('taking the default parameters')
	Nx=25; 				# size along x 
	Ny=25; 				# size along y 
	radius=8;			# radius of central lead 
	Niter=1500; 		# number of iterations to perform
phi=zeros([Ny,Nx])               # phi is potential array
x=linspace(-(Nx-1)/2,(Nx-1)/2,Nx)
y=linspace((Ny-1)/2,-(Ny-1)/2,Ny) 
X,Y=meshgrid(x,y)
xy=where(X*X+Y*Y<(radius*radius))
phi[xy]=1						# Making potetial at points where x2+y2<c2
electrode=phi.copy()
def update_phi(phi):			# updating potential array i.e., average of neighbour potentials
    phi[1:-1,1:-1]=0.25*(phi[1:-1,0:-2]+phi[1:-1,2:]+phi[0:-2,1:-1]+phi[2:,1:-1])
    return phi
def boundaries_phi(phi):		# Boundary Conditions
    phi[:,0]=phi[:,1]
    phi[:,-1]=phi[:,-2] 
    phi[0,:]=phi[1,:] 
    phi[-1,:]=0
    return phi
errors=zeros(Niter)				# Error array for storing amount of change for each iteration
# Getting final potential array after Niter iterations at corresponding points
for k in range(Niter):
    oldphi=phi.copy()
    update_phi(phi)
    boundaries_phi(phi)
    phi[xy]=1
    errors[k]=((abs(phi-oldphi)).max())			# Error array for storing amount of change for each iteration
A=[]
for i in range(Niter):
    A.append([1,i])
p1=lstsq(c_[A],c_[log(errors)],rcond=1)[0]				# Getting bestfit of A and B for y=Aexp(Bx) from the start of iteration
p2=lstsq(c_[A[500:]],c_[log(errors[500:])],rcond=1)[0]	# Getting bestfit of A and B for y=Aexp(Bx) from the start of 500th iteration
errfit1=e**dot(c_[A],p1)							# Storing errors got by putting bestfit A and B values of y=Aexp(Bx) from the start of iteration
errfit2=e**dot(c_[A],p2)							# Storing errors got by putting bestfit A and B values of y=Aexp(Bx) from the start of 500th iteration
jx=zeros([Ny,Nx])
jy=zeros([Ny,Nx])
jx[1:-1,1:-1]=0.5*(phi[1:-1,0:-2]-phi[1:-1,2:])		# Getting values of current density in x direction by taking difference in potential in x
jy[1:-1,1:-1]=0.5*(phi[0:-2,1:-1]-phi[2:,1:-1])		# Getting values of current density in y direction by taking difference in potential in y
figure(1)
semilogy(c_[A][0:,1:2][::50],errors[::50],'r.',markersize=14,label='errors')
semilogy(c_[A][0:,1:2][::50],errfit1[::50],'g.',markersize=11,label='fit1')
semilogy(c_[A][0:,1:2][::50],errfit2[::50],'b.',markersize=8,label='fit2')			# Plotting of errors in semilogy plot
legend(loc='upper right')
xlabel(r'$n$',size=10)
ylabel(r'$errors$',size=10)
figure(6)
loglog(c_[A][0:,1:2][::50],errors[::50],'g.',label='errorslog')
legend(loc='upper right')
xlabel(r'$n$',size=10)
ylabel(r'$errors$',size=10)
figure(2)  
ax=p3.Axes3D(figure(2)) # Axes3D is the means to do a surface plot 
title('The 3-D surface plot of the potential') 						
surf = ax.plot_surface(X,Y, phi, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)				# Surface plot
colorbar(surf)
xlabel(r'$x$',size=10)
ylabel(r'$y$',size=10)
figure(3)
con=contourf(X,Y,phi,Nx)					# Plotting contour of potential phi
plot(X[xy],Y[xy],'r.')
colorbar(con)
figure(4)
con2=contourf(X,Y,electrode)
colorbar(con2)					# plotting 1 volt potential
figure(5)
title('The vector plot of the current flow') 
plot(X[xy],Y[xy],'r.')
quiver(X,Y,-jx,-jy,headlength=3,headwidth=3,scale=4)						# plotting direction of currents from higher potential to lower potential
show()