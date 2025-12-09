import numpy as np
import matplotlib.pyplot as plt

#load file
data = np.loadtxt('../outputs/Legendre/legendre_polynomials.txt', delimiter=',')

#plot data
x = data[:,0]
P0 = data[:,1]
P1 = data[:,3]
P2 = data[:,5]
P3 = data[:,7]
P4 = data[:,9]
P5 = data[:,11]

plt.plot(x, P0, label='P0')
plt.plot(x, P1, label='P1')
plt.plot(x, P2, label='P2')
plt.plot(x, P3, label='P3')
plt.plot(x, P4, label='P4')
plt.plot(x, P5, label='P5')
plt.xlabel('x')
plt.ylabel('P_n(x)')
plt.title('Legendre Polynomials')
plt.legend()
plt.grid()
plt.savefig('../outputs/Legendre/legendre_polynomials.png')
plt.clf()

dP0 = data[:,2]
dP1 = data[:,4]
dP2 = data[:,6]
dP3 = data[:,8]
dP4 = data[:,10]
dP5 = data[:,12]
plt.plot(x, dP0, label="P0'")
plt.plot(x, dP1, label="P1'")
plt.plot(x, dP2, label="P2'")
plt.plot(x, dP3, label="P3'")
plt.plot(x, dP4, label="P4'")
plt.plot(x, dP5, label="P5'")
plt.xlabel('x')
plt.ylabel("P_n'(x)")
plt.title("Derivatives of Legendre Polynomials")
plt.legend()
plt.grid()
plt.savefig('../outputs/Legendre/legendre_polynomials_derivatives.png')
plt.clf()
