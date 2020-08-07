#																		    	  			END-SEMESTER EXAMINATION
#																			        		 NITHIN BABU [EE18B021]

# Importing the required modules.	
from pylab import *

# Declaring the function to solve the Laplace's Equation and calculate the Potential Matrix.
def potential_calculator(M, N, delta, k, accuracy, N_max, graph = 'No'):

# Declaring all the necessary variables and making the meshgrid.
	E_r = 2; N_iter = 0
	Ly = arange(0, delta*M, delta)
	Lx = arange(0, delta*N, delta)

	X,Y = meshgrid(Ly,Lx)

	error = zeros(N_max)
	phi = zeros((M,N))

# This for loop is run for N_max number of times so that the potential matrix converges.
	for p in range(N_max):

		N_iter = N_iter + 1

		oldphi = phi.copy()

# This is the parallelize algorithm followed. The matrix is first divided into two parts and calculated accordingly. 
		phi[1:k,1:-1] = 0.25*(phi[1:k,0:-2] + phi[1:k,2:] + phi[0:k-1,1:-1] + phi[2:k+1,1:-1])
		phi[k+1:-1,1:-1] = 0.25*(phi[k+1:-1,0:-2] + phi[k+1:-1,2:] + phi[k:-2,1:-1] + phi[k+2:,1:-1])

# The boundary conditions are also being assigned here.
		phi[k,1:-1] = (E_r*phi[k-1,1:-1] + phi[k+1,1:-1])/(1 + E_r)
		phi[-1,1:-1] = 1.0

# The error matrix is updated so that, once the error is less than the accuracy, the loop is broken.
		error[p] = abs(oldphi - phi).max() 

		if error[p] <= accuracy:
			break;

# A spacial contour plot for the potential matrix.
	if graph == 'Yes':

		figure(0)	
		contourf(Y,X,phi.T)	
		colorbar()
		title("The Potential Diagram")
		xlabel(r'$X(cm)\rightarrow$')
		ylabel(r'$Y(cm)\rightarrow$')
		grid(True)

# The function returns the potential matrix, the error vector and the number of iterations necessary.
	return phi, error, N_iter


# Declaring the function to calculate the two components of the Electric field from the potential Matrix.
def field_calculator(phi, delta):

# The below piece of code will calculate the components of Electric field at the centre of the mesh cells.
	Ex = 100*(phi[1:,0:-1] - phi[1:,1:])/(delta)
	Ey = 100*(phi[0:-1,1:] - phi[1:,1:])/(delta)

	return Ex, Ey	

# Declaring the function to calculate the top plate charge and the fluid charge using the Electric field components.
def charge_calculator(Ex,Ey,delta,Lz):

# The below piece of code is an approximate to the integration required to get the charge from Gauss Law.
	Qtop = abs(sum(Epsilon*Ey[-1,:]*delta))*Lz*1e-4

	Qfluid = -abs(sum(Epsilon*Ey[0,:]*delta))*Lz*1e-4 - (abs(sum(Epsilon*Ex[0:k-1,-1]*delta))*Lz*1e-4) 

	return Qtop, Qfluid

# Declaring all the necessary variables required for answering the other parts of the question.
# Let us assume the other dimension Lz to be 10cm.
Lz = 10
Epsilon = 8.85e-12
M = 21; N = 11; delta = 1
k = 10; accuracy = 1e-5; N_max = 1000

# Calling the function potential_calculator() to find the potential matrix, error vector and the number of iterations carried out.
phi, error, N_iter = potential_calculator(M,N,delta,k,accuracy,N_max,'Yes')

print("The number of iterations carried out is: ",N_iter)

# Construction of the meshgrid.
Lx = arange(0.5, delta*(N-1), delta)
Ly = arange(0.5, delta*(M-1), delta)
X,Y = meshgrid(Ly,Lx)

#Calculating the values of the two charges for different values of h.
Qtop = []
Qfluid = []
for h in arange(0.1,1,0.1):

# Calling all the functions one by one to calculate the 2 charge values.
	k = int(h*(M-1))
	phi,_,_ = potential_calculator(M,N,delta,k,accuracy,N_max)
	Ex, Ey = field_calculator(phi, delta)
	Qt, Qf = charge_calculator(Ex, Ey, delta, Lz)
	Qtop.append(Qt)
	Qfluid.append(Qf)

# The below piece of code is to show that Dn remains continuous before and after the fluid boundary.
k = 10
phi,_,_ = potential_calculator(M,N,delta,k,accuracy,N_max)
Ex, Ey = field_calculator(phi, delta)

# Calculating the values of Dn0 (Dn after the boundary), and Dn1 (Dn before the boundary)
Dn0 = abs(Ey[k-1,:-1]*2*Epsilon); Dn1 = abs(Ey[k,:-1]*Epsilon)
err_D = abs((Dn1 - Dn0))/Dn0;

# Printing the two values of Dn along with the relative mean error.
print("\nThe average values of Dn0 and Dn1 at the surface Yk are: ")
print("\nDno: ",mean(Dn0))
print("\nDn1: ",mean(Dn1))
print("\nThe relative mean error is :",mean(err_D))  

if mean(err_D)<=1e-15:
	print("Hence, Dn is continuous.")

else:
	print("Hence, Dn is not continuous.")	


# The below piece of code is to check if Snell's Law is valid at the boundary.

# The Snell's terms are calculated using Ex, Ey.
tetha_1 = arctan(Ey[k-1][5]/Ex[k-1][5]); tetha_2 = arctan(Ey[k][5]/Ex[k][5])
Snell_term_1 = sqrt(2)*sin(tetha_1); Snell_term_2 = sin(tetha_2)

err_snell = abs((Snell_term_2 - Snell_term_1)/Snell_term_1)

print("\nThe change in the angle of the Electric field is :",abs(tetha_2-tetha_1))

# Printing the relative error in between the Snell terms.
print("\nThe Snell's Law relative error is: ",err_snell)

if err_snell<=1e-15:
	print("Hence, Snell's Law is valid.")

else:
	print("Hence, Snell's Law is invalid.")	

# Plot of the Error vector vs the iteration number.
figure(1)
plot(arange(N_iter),error[where(error!=0)],label='Error')
title("The Error Plot")
xlabel(r'$N_{max}\rightarrow$')
ylabel(r'$Error\rightarrow$')
legend()
grid(True)

# The Quiver plot of the spacial Electric field.
figure(2)
quiver(Y,X,Ex.T,Ey.T)
title("The Electric Field Plot")
xlabel(r'$X(cm)\rightarrow$')
ylabel(r'$Y(cm)\rightarrow$')
grid(True)

# The plot of the top plate charge vs h/Ly.
figure(3)
plot(arange(0.1,1,0.1),Qtop,label='Top Plate Charge')
title("The Top plate Charge Plot")
xlabel(r'$h/L_y\rightarrow$')
ylabel(r'$Charge\rightarrow$')
legend()
grid(True)

# The plot of the fluid charge vs h/Ly.
figure(4)
plot(arange(0.1,1,0.1),Qfluid,label='Fluid Charge')
title("The Fluid Charge Plot")
xlabel(r'$h/L_y\rightarrow$')
ylabel(r'$Charge\rightarrow$')
legend()
grid(True)

show()
