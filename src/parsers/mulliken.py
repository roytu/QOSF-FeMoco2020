
import numpy as np
import matplotlib.pyplot as plt

def read_mulliken_file(fname):
    with open(fname, "r") as f:
        lines = f.readlines()

    mulliken = []
    for line in lines:
        elems = line.split()
        atom_num = elems[0]
        atom_type = elems[1]
        atom_orb = elems[2]
        atom_pop = float(elems[3])

        orbname = f"{atom_num} {atom_type} {atom_orb}"

        mulliken.append((orbname, atom_pop))

    # Sort by population
    mulliken = sorted(mulliken, key=lambda xs: xs[1], reverse=True)

    # Remove orbitals greater than threshold
    THRESHOLD = 1.995
    mulliken = list(filter(lambda xs: xs[1] <= THRESHOLD, mulliken))

    # Pretty print
    for (orbname, atom_pop) in mulliken:
        print(f"{orbname}\t\t{atom_pop}")
    print(f"Number of active orbitals: {len(mulliken)}")

    # Graph distribution
    #pops = [pop for (_, pop) in mulliken]

    #plt.title("Mulliken occupation number distribution")
    #plt.hist(pops, bins=1000)
    #plt.show()

if __name__ == "__main__":
    read_mulliken_file("mulliken.txt")
