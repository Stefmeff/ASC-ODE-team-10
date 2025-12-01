import numpy as np

data = []


for steps in [50, 100, 150, 200]:
    

    filename1 = f'../outputs/MassSpring/ExplicitEuler/steps{steps}.txt'
    filename2 = f'../outputs/MassSpring/ImprovedEuler/steps{steps}.txt'
    filename3 = f'../outputs/MassSpring/ImplicitEuler/steps{steps}.txt'
    filename4 = f'../outputs/MassSpring/CrankNicolson/steps{steps}.txt'
    filename5 = f'../outputs/MassSpring/RungeKutta/steps{steps}.txt'
    #  filename6 = f'../outputs/ImplicitRungeKutta/steps{steps}.txt'
    d1 = np.loadtxt(filename1, usecols=(0, 1, 2))
    d2 = np.loadtxt(filename2, usecols=(0, 1, 2))
    d3 = np.loadtxt(filename3, usecols=(0, 1, 2))
    d4 = np.loadtxt(filename4, usecols=(0, 1, 2))
    d5 = np.loadtxt(filename5, usecols=(0, 1, 2))
    # d6 = np.loadtxt(filename6, usecols=(0, 1, 2))
    data = [d1, d2, d3, d4, d5]


    import matplotlib.pyplot as plt


    plt.plot(data[0][:,0], data[0][:,1], label='position-explicit Euler')
    plt.plot(data[1][:,0], data[1][:,1], label='position-improved Euler')
    plt.plot(data[2][:,0], data[2][:,1], label='position-implicit Euler')
    plt.plot(data[3][:,0], data[3][:,1], label='position-Crank-Nicolson')
    plt.xlabel('time')
    plt.ylabel('value')
    plt.title('Mass-Spring System Time Evolution')
    plt.legend()
    plt.grid()
    plt.savefig(f"../outputs/plots/massspring_position_steps{steps}.png")
    plt.clf()

    # Plot Runge Kutta methods
    plt.plot(data[4][:,0], data[4][:,1], label='position-Explicit Runge Kutta')
    # plt.plot(data[5][:,0], data[5][:,1], label='position-Implicit Runge Kutta')
    plt.xlabel('time')
    plt.ylabel('value')
    plt.title('Mass-Spring System Time Evolution')
    plt.legend()
    plt.grid()
    plt.savefig(f"../outputs/plots/Runge_massspring_position_steps{steps}.png")
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
    plt.savefig(f"../outputs/plots/massspring_phase_steps{steps}.png")
    plt.clf()

    # Plot Runge Kutta methods phase plot
    plt.plot(data[4][:,1], data[4][:,2], label='phase plot-Explicit Runge Kutta')
    # plt.plot(data[5][:,1], data[5][:,2], label='phase plot-Implicit Runge Kutta')
    plt.xlabel('position')
    plt.ylabel('velocity') 
    plt.title('Mass-Spring System Phase Plot')
    plt.legend()
    plt.grid()
    plt.savefig(f"../outputs/plots/Runge_massspring_phase_steps{steps}.png")
    plt.clf()