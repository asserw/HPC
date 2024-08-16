import matplotlib.pyplot as plt

def plot_energy(filenames, timesteps):
    plt.figure()

    for filename, timestep in zip(filenames, timesteps):
        times = []
        energies = []
        with open(filename) as f:
            for step, line in enumerate(f):
                energies.append(float(line.strip()))
                times.append(step * timestep)
        plt.plot(times, energies, label=f'Time step: {timestep}')

    plt.xlabel('Time')
    plt.ylabel('Total Energy')
    plt.title('Total Energy vs Time for Different Time Steps')
    plt.legend()
    plt.savefig("total_energy_combined_plot.png")

filenames = ["energy_0.001000.txt", "energy_0.004000.txt", "energy_0.008000.txt", "energy_0.010000.txt"]
timesteps = [0.001, 0.004, 0.008, 0.010]

plot_energy(filenames, timesteps)
