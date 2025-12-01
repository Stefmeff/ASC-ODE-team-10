import numpy as np

data = []


for steps in [50, 100, 150, 200]:
    

    filename1 = f'../outputs/ElectricNetwork/ExplicitEuler/steps{steps}.txt'
    filename2 = f'../outputs/ElectricNetwork/ImplicitEuler/steps{steps}.txt'
    filename3 = f'../outputs/ElectricNetwork/CrankNicolson/steps{steps}.txt'    
    d1 = np.loadtxt(filename1, usecols=(0, 1, 2))
    d2 = np.loadtxt(filename2, usecols=(0, 1, 2))
    d3 = np.loadtxt(filename3, usecols=(0, 1, 2))
    data = [d1, d2, d3]


    import matplotlib.pyplot as plt


    plt.plot(data[0][:,0], data[0][:,1], label='position-explicit Euler')
    plt.plot(data[1][:,0], data[1][:,1], label='position-implicit Euler')
    plt.plot(data[2][:,0], data[2][:,1], label='position-Crank-Nicolson')
    plt.xlabel('time')
    plt.ylabel('voltage U_c')
    plt.title('Electric Network Time Evolution')
    plt.legend()
    plt.grid()
    plt.savefig(f"../outputs/ElectricNetwork/plots/voltage_plot_steps{steps}.png")
    plt.clf()

