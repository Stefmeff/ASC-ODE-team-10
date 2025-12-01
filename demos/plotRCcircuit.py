import numpy as np

data = []


steps = 100


filename1 = f'../outputs/ElectricNetwork/ExplicitEuler/steps{steps}.txt'
filename2 = f'../outputs/ElectricNetwork/ImplicitEuler/steps{steps}.txt'
filename3 = f'../outputs/ElectricNetwork/CrankNicolson/steps{steps}.txt'   
d1 = np.loadtxt(filename1, usecols=(0, 1, 2))
d2 = np.loadtxt(filename2, usecols=(0, 1, 2))
d3 = np.loadtxt(filename3, usecols=(0, 1, 2))
data = [d1, d2, d3]

import matplotlib.pyplot as plt


plt.plot(data[0][:,0], data[0][:,1], label='explicit Euler')
plt.plot(data[1][:,0], data[1][:,1], label='implicit Euler')
plt.plot(data[2][:,0], data[2][:,1], label='Crank-Nicolson')
plt.xlabel('time')
plt.ylabel('voltage U_c')
plt.title('Electric Network Time Evolution')
plt.legend()
plt.grid()
plt.savefig(f"../outputs/ElectricNetwork/plots/voltage_plot_steps{steps}.png")
plt.clf()

# Implicit Runge Kutta method
filename4 = f'../outputs/ElectricNetwork/ImplicitRungeKutta/stages:2steps:{steps}.txt'
filename5 = f'../outputs/ElectricNetwork/ImplicitRungeKutta/stages:3steps:{steps}.txt'
filename6 = f'../outputs/ElectricNetwork/ImplicitRungeKutta/stages:5steps:{steps}.txt'
d4 = np.loadtxt(filename4, usecols=(0, 1, 2))
d5 = np.loadtxt(filename5, usecols=(0, 1, 2))
d6 = np.loadtxt(filename6, usecols=(0, 1, 2))
data_rk = [d4, d5, d6]
plt.plot(data_rk[0][:,0], data_rk[0][:,1], label='2 stages')
plt.plot(data_rk[1][:,0], data_rk[1][:,1], label='3 stages')
plt.plot(data_rk[2][:,0], data_rk[2][:,1], label='5 stages')
plt.xlabel('time')
plt.ylabel('voltage U_c')
plt.title('Electric Network Implicit Runge Kutta Time Evolution')
plt.legend()
plt.grid()
plt.savefig(f"../outputs/ElectricNetwork/plots/ImplicitRunge_voltage_plot_steps{steps}.png")
plt.clf()

# Explicit Runge Kutta method
filename7 = f'../outputs/ElectricNetwork/ExplicitRungeKutta/stages:2steps:{steps}.txt'
filename8 = f'../outputs/ElectricNetwork/ExplicitRungeKutta/stages:3steps:{steps}.txt'
filename9 = f'../outputs/ElectricNetwork/ExplicitRungeKutta/stages:5steps:{steps}.txt'
d7 = np.loadtxt(filename7, usecols=(0, 1, 2))
d8 = np.loadtxt(filename8, usecols=(0, 1, 2))
d9 = np.loadtxt(filename9, usecols=(0, 1, 2))
data_rk_exp = [d7, d8, d9]
plt.plot(data_rk_exp[0][:,0], data_rk_exp[0][:,1], label='2 stages')
plt.plot(data_rk_exp[1][:,0], data_rk_exp[1][:,1], label='3 stages')
plt.plot(data_rk_exp[2][:,0], data_rk_exp[2][:,1], label='5 stages')
plt.xlabel('time')
plt.ylabel('voltage U_c')
plt.title('Electric Network Explicit Runge Kutta Time Evolution')
plt.legend()
plt.grid()
plt.savefig(f"../outputs/ElectricNetwork/plots/ExplicitRunge_voltage_plot_steps{steps}.png")
plt.clf()   

