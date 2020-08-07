#																		ASSIGNMENT 3: FITTING DATA TO MODELS														
#																				NITHIN BABU [EE18B021]			

# Importing the required modules.
from pylab import *
import scipy.special as sp

# Defining the actual function. 
def g(t,A,B):
	return A*sp.jn(2,t) + B*t

# Loading data from "fitting.dat"
data = loadtxt("fitting.dat")

# Extracting time and the first set of data.
t = data[:,0]
c = data[:,1]

# Calculating the actual set of functional values.
y0 = g(t,1.05,-0.105)	
stdev0 = 0.10
stdev = logspace(-1,-3,9)
n,m = data.shape

# Calculating the M matrix asked in Q.6.
M = empty((n,2))
for i in range(n):
	M[i] = (sp.jn(2,t[i]),t[i])

A0 = array([1.05,-0.105])
y = dot(M,A0)

if array_equal(y,g(t,1.05,-0.105)):
	print("Both solutions for Q.6 are equal.")
else:
	print("Both solutions for Q.6 are not equal.")	

# The below lines of code is used to calculate the mean squared error for different values of A and B.
A = array([0.1*i for i in range(21)])
B = array([-0.2+0.01*i for i in range(21)])
E = zeros((21,21))
for i in range(21):
	for j in range(21):
		for k in range(n):
			E[i][j] += ((c[k]-g(t[k],A[i],B[j]))**2)/n

# This is used to plot the contour. 
X,Y = meshgrid(A,B)

# The below part of code is used to calculate the error in the estimate of A and B.
Ea = empty((9,1))
Eb = empty((9,1))

for j in range(9):
	
	AB = linalg.lstsq(M,data[:,j+1],rcond=None)
	Ea[j] = abs(AB[0][0]-A0[0])
	Eb[j] = abs(AB[0][1]-A0[1])

# This will print the actual plot along with nine other noisy plots.
figure(0)
for i in range(1,10):
	plot(t,data[:,i],label="σ=%.3f"%stdev[i-1])
plot(t,y0,label="True Value",color='black',linewidth=3)	
title("Figure 0",size=20)	
xlabel(r'$t\rightarrow$',size=20)
ylabel(r'$f(t)+n\rightarrow$',size=20)
grid(True)
legend()

# This will plot the error plot in the first set of data as compared to the actual plot.
figure(1)
plot(t,y0,label="True Value",color='black',linewidth=3)
errorbar(t[::5],c[::5],stdev0,fmt='ro',label='Noise')
title("Data with Error Bars",size=20)	
xlabel(r'$t\rightarrow$',size=20)
ylabel(r'$f(t)+n\rightarrow$',size=20)
grid(True)
legend()

# This will plot the contour of mean squared error for different values of A and B.
figure(2)
Contour=contour(X,Y,E,25)
clabel(Contour,[0.02,0.04,0.06,0.08],inline=1)
title("Contour Plot",size=20)	
xlabel(r'$A\rightarrow$',size=20)
ylabel(r'$B\rightarrow$',size=20)
grid(True)

# This will plot the variation of error in the estimate of A and B with respect to noise.
figure(3)
plot(stdev,Ea,label='Aerr',marker='o',linestyle='dashed')
plot(stdev,Eb,label='Berr',marker='o',linestyle='dashed')
title("Variation of Error with Noise",size=20)	
xlabel('Noise standard deviation',size=20)
ylabel('MS error',size=20)
grid(True)
legend()

# This will plot the variation of error in the estimate of A and B with respect to noise in log scale.
figure(4)
loglog(stdev,Ea,'ro',label='Aerr',)
errorbar(logspace(-1, -3, 9), Ea, std(Ea), fmt='ro')
loglog(stdev,Eb,'go',label='Berr')
errorbar(logspace(-1, -3, 9), Eb, std(Eb), fmt='go')
title("Variation of Error with Noise",size=20)	
xlabel('Noise standard deviation',size=20)
ylabel('MS error',size=20)
grid(True)
legend()

# This command will display all the graphs defined above.
show()
