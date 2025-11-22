import numpy as np

data = []


for steps in [50, 100, 150, 200]:
    

    filename1 = f'../outputs/ExplicitEuler/steps{steps}.txt'
    filename2 = f'../outputs/ImprovedEuler/steps{steps}.txt'
    filename3 = f'../outputs/ImplicitEuler/steps{steps}.txt'
    filename4 = f'../outputs/CrankNicolson/steps{steps}.txt'
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
    plt.ylabel('value')
    plt.title('Mass-Spring System Time Evolution')
    plt.legend()
    plt.grid()
    plt.savefig(f"../outputs/plots/massspring_position_steps{steps}.png")
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

