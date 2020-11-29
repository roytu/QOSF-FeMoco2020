
import re
import numpy as np
import matplotlib.pyplot as plt


def graph():
    with open("rohf_progress.out", "r") as f:
        lines = f.readlines()
    
    cycles = []
    energies = []
    for line in lines:
        m = re.search("cycle=\s([0-9]*)\sE=\s([0-9\.-]*)\s", line)
        cycle = int(m.group(1))
        energy = float(m.group(2))

        cycles.append(cycle)
        energies.append(energy)

    cycles = np.array(cycles)
    energies = np.array(energies)

    #np.savetxt("rohf_energies.csv", [cycles, energies])
    
    plt.figure()
    plt.title("ROHF Energies (FeMoco)")
    plt.xlabel("Cycles")
    plt.ylabel("Energies")
    plt.plot(cycles, energies)
    plt.savefig("rohf_energies.png")
    plt.show()

if __name__ == "__main__":
    graph()
