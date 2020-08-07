#																 ASSIGNMENT 7: CIRCUIT ANALYSIS USING SYMPY AND LAPLACE TRANSFORMS
#																			           NITHIN BABU [EE18B021]
# Importing the necessary modules.
from sympy import *
import scipy.signal as sp
import pylab as p

# Declaring the sympy function for a lowpass filter.
def lowpass(R1,R2,C1,C2,G,Vi):
	s = symbols('s')
# Creating the matrices and solving them to get the output voltage.	
	A = Matrix([[0,0,1,-1/G],
		[-1/(1+s*R2*C2),1,0,0],
		[0,-G,G,1],
		[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
	b = Matrix([0,0,0,-Vi/R1])
	V = A.inv()*b
	return A,b,V

# Declaring the sympy function for a highpass filter.
def highpass(R1,R3,C1,C2,G,Vi):
    s = symbols("s")
# Creating the matrices and solving them to get the output voltage.	
    A = Matrix([[0,-1,0,1/G],
        [s*C2*R3/(s*C2*R3+1),0,-1,0],
        [0,G,-G,1],
        [-s*C2-1/R1-s*C1,0,s*C2,1/R1]])
    b = Matrix([0,0,0,-Vi*s*C1])
    V = A.inv()*b
    return A,b,V

# Creating a function to convert a sympy function into a version that is understood by sp.signal  	 
def sympyToTrFn(Y):
# The below line will simplify the expression and give it in Nr/Dr form.	
    Y = expand(simplify(Y))
# The following lines will give the proper form for numerator and denomerator.    
    n,d = fraction(Y)
    n,d = Poly(n,s), Poly(d,s)
    num,den = n.all_coeffs(), d.all_coeffs()
    num,den = [float(f) for f in num], [float(f) for f in den]
# This will calculate the function that is understood by sp.signal    
    H = sp.lti(num,den)
    return H

# Reassigning the value of pi.
pi = p.pi

# The below piece of code will calculate the transfer function for the lowpass filter.
s = symbols('s')
A,b,V = lowpass(10000,10000,1e-9,1e-9,1.586,1)
Vo = V[3]
H = sympyToTrFn(Vo)
ww = p.logspace(0,8,801)
ss = 1j*ww
hf = lambdify(s,Vo,'numpy')
v = hf(ss)

# These lines of code will calculate the step response for the lowpass filter.
A1,b1,V1 = lowpass(10000,10000,1e-9,1e-9,1.586,1/s)
Vo1 = V1[3]
H1 = sympyToTrFn(Vo1)
t,y1 = sp.impulse(H1,None,p.linspace(0,5e-3,10000))

# The response is also calculated for sum of sinusoids.
vi = p.sin(2000*pi*t) + p.cos(2e6*pi*t)
t,y2,svec = sp.lsim(H,vi,t)

# The below piece of code will calculate the transfer function for the highpass filter.
A3,b3,V3 = highpass(10000,10000,1e-9,1e-9,1.586,1)
Vo3 = V3[3]
H3 = sympyToTrFn(Vo)
hf3 = lambdify(s,Vo3,'numpy')
v3 = hf3(ss)

# The output response is calculated when the input is a damped sinusoid for both high and low frequency.
# High frequency.
damping_factor = -500
t2 = p.linspace(0,1e-2,1e5) 
vi4_1 = p.exp(damping_factor*t2)*p.cos(2e6*pi*t2)
t2,y4_1,svec = sp.lsim(H3,vi4_1,t2)

# Low frequency.
t3 = p.linspace(0,1e-2,1e5)
vi4_2 = p.exp(damping_factor*t3)*p.cos(2e3*pi*t3)
t3,y4_2,svec = sp.lsim(H3,vi4_2,t3)

# The step response is calculated for a highpass filter.
A5,b5,V5 = highpass(10000,10000,1e-9,1e-9,1.586,1/s)
Vo5 = V5[3]
H5 = sympyToTrFn(Vo5)
t5,y5 = sp.impulse(H5,None,p.linspace(0,5e-3,10000))

# The plot for Magnitude of transfer function of a lowpass filter.
p.figure(0)
p.loglog(ww,abs(v),lw=2)
p.title(r"$|H(j\omega)|$ for lowpass filter")
p.xlabel(r'$\omega\rightarrow$')
p.ylabel(r'$|H(j\omega)|\rightarrow$')
p.grid(True)

# The plot for step response of a lowpass filter.
p.figure(1)
p.plot(t,y1)
p.title(r"Step Response for low pass filter")
p.xlabel(r'$t\rightarrow$')
p.ylabel(r'$V_o(t)\rightarrow$')
p.grid(True)

# The plot for output response for sum of sinusoids of a lowpass filter.
p.figure(2)
p.plot(t,y2)
p.title(r"Output voltage for sum of sinusoids")
p.xlabel(r'$t\rightarrow$')
p.ylabel(r'$V_o(t)\rightarrow$')
p.grid(True)

# The plot for Magnitude of transfer function of a highpass filter.
p.figure(3)
p.loglog(ww,abs(v3),lw=2)
p.title(r"$|H(j\omega)|$ for highpass filter")
p.xlabel(r'$\omega\rightarrow$')
p.ylabel(r'$|H(j\omega)|\rightarrow$')
p.grid(True)

# The plot for high frequency damped sinusoids.
p.figure(4)
p.plot(t2,vi4_1)
p.title(r"High frequency damped sinusoid ")
p.xlabel(r'$t\rightarrow$')
p.ylabel(r'$V_o(t)\rightarrow$')
p.grid(True)

# The plot for low frequency damped sinusoids.
p.figure(5)
p.plot(t3,vi4_2)
p.title(r"Low frequency damped sinusoid")
p.xlabel(r'$t\rightarrow$')
p.ylabel(r'$V_o(t)\rightarrow$')
p.grid(True)

# The plot for high frequency damped sinusoid response from highpass filter.
p.figure(6)
p.plot(t2,y4_1)
p.title(r"High frequency damped sinusoid response from High Pass filter")
p.xlabel(r'$t\rightarrow$')
p.ylabel(r'$V_o(t)\rightarrow$')
p.grid(True)

# The plot for low frequency damped sinusoid response from highpass filter.
p.figure(7)
p.plot(t3,y4_2)
p.title(r"Low frequency damped sinusoid response from High Pass filter")
p.xlabel(r'$t\rightarrow$')
p.ylabel(r'$V_o(t)\rightarrow$')
p.grid(True)

# The plot for step response of a highpass filter.
p.figure(8)
p.plot(t5,y5)
p.title(r"Step Response for high pass filter")
p.xlabel(r'$t\rightarrow$')
p.ylabel(r'$V_o(t)\rightarrow$')
p.grid(True)

p.show()	





