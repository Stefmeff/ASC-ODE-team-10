import numpy as np

data = []


for steps in [50, 100, 150, 200]:
    

    filename1 = f'../outputs/MassSpring/ExplicitEuler/steps{steps}.txt'
    filename2 = f'../outputs/MassSpring/ImprovedEuler/steps{steps}.txt'
    filename3 = f'../outputs/MassSpring/ImplicitEuler/steps{steps}.txt'
    filename4 = f'../outputs/MassSpring/CrankNicolson/steps{steps}.txt'
    d1 = np.loadtxt(filename1, usecols=(0, 1, 2))
    d2 = np.loadtxt(filename2, usecols=(0, 1, 2))
    d3 = np.loadtxt(filename3, usecols=(0, 1, 2))
    d4 = np.loadtxt(filename4, usecols=(0, 1, 2))
    data = [d1, d2, d3, d4]


    import matplotlib.pyplot as plt


    plt.plot(data[0][:,0], data[0][:,1], label='position-explicit Euler')
    plt.plot(data[1][:,0], data[1][:,1], label='position-improved Euler')
    plt.plot(data[2][:,0], data[2][:,1], label='position-implicit Euler')
    plt.plot(data[3][:,0], data[3][:,1], label='position-Crank-Nicolson')
    plt.xlabel('time')
    plt.ylabel('position')
    plt.title('Mass-Spring System Time Evolution')
    plt.legend()
    plt.grid()
    plt.savefig(f"../outputs/MassSpring/plots/massspring_position_steps{steps}.png")
    plt.clf()





    plt.plot(data[0][:,1], data[0][:,2], label='phase plot-explicit Euler')
    plt.plot(data[1][:,1], data[1][:,2], label='phase plot-improved Euler')
    plt.plot(data[2][:,1], data[2][:,2], label='phase plot-implicit Euler')
    plt.plot(data[3][:,1], data[3][:,2], label='phase plot-Crank-Nicolson')
    plt.xlabel('position')
    plt.ylabel('velocity')
    plt.title('Mass-Spring System Phase Plot')
    plt.legend()
    plt.grid()
    plt.savefig(f"../outputs/MassSpring/plots/massspring_phase_steps{steps}.png")
    plt.clf()

# Plot Implicit Runge Kutta methods
filename5 = f'../outputs/MassSpring/ImplicitRungeKutta/stages:2steps:150.txt'
filename6 = f'../outputs/MassSpring/ImplicitRungeKutta/stages:3steps:150.txt'
filename7 = f'../outputs/MassSpring/ImplicitRungeKutta/stages:5steps:150.txt'
data5 = np.loadtxt(filename5, usecols=(0, 1, 2))
data6 = np.loadtxt(filename6, usecols=(0, 1, 2))
data7 = np.loadtxt(filename7, usecols=(0, 1, 2))
data = [data5, data6, data7]    

plt.plot(data[0][:,0], data[0][:,1], label='2 stages')
plt.plot(data[1][:,0], data[1][:,1], label='3 stages')
plt.plot(data[2][:,0], data[2][:,1], label='5 stages')
plt.xlabel('time')
plt.ylabel('position')
plt.title('Mass-Spring System Implicit Runge Kutta Time Evolution')
plt.legend()
plt.grid()
plt.savefig(f"../outputs/MassSpring/plots/ImpRunge_massspring_position_steps150.png")
plt.clf()

#Plot Explicit Runge Kutta methods
filename8 = f'../outputs/MassSpring/ExplicitRungeKutta/stages:2steps:150.txt'
filename9 = f'../outputs/MassSpring/ExplicitRungeKutta/stages:3steps:150.txt'
filename10 = f'../outputs/MassSpring/ExplicitRungeKutta/stages:5steps:150.txt'
data8 = np.loadtxt(filename8, usecols=(0, 1, 2))
data9 = np.loadtxt(filename9, usecols=(0, 1, 2))
data10 = np.loadtxt(filename10, usecols=(0, 1, 2))
data = [data8, data9, data10]           
plt.plot(data[0][:,0], data[0][:,1], label='2 stages')
plt.plot(data[1][:,0], data[1][:,1], label='3 stages')
plt.plot(data[2][:,0], data[2][:,1], label='5 stages')  
plt.xlabel('time')
plt.ylabel('position')
plt.title('Mass-Spring System Explicit Runge Kutta Time Evolution')
plt.legend()
plt.grid()
plt.savefig(f"../outputs/MassSpring/plots/ExpRunge_massspring_position_steps150.png")
plt.clf()