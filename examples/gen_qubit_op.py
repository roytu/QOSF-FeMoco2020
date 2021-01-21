
import sys
import os
import pathlib
import time

import logging
logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)
#logging.getLogger().setLevel(logging.INFO)

import numpy as np

# ----- Import modified qiskit -----

sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), "lib")))

from qiskit.chemistry.drivers import PySCFDriver, UnitsType
from qiskit.chemistry import QMolecule, FermionicOperator
from qiskit.chemistry.core import Hamiltonian, TransformationType, QubitMappingType
from qiskit.chemistry.transformations import FermionicTransformation

class GenQubitOp(object):
    def __init__(self):
        # -- Configure this --
        self.filename = "hdf5_files/lih_sto3g_(0,0).hdf5"
        self.geometry_filename = "molecules/ICS.xyz"
        self.charge = 0
        self.spin = 0
        self.freeze_list = [2]
        self.remove_list = [1, 9]
        self.map_type = "jordan_wigner" # parity, jordan_wigner, or bravyi_kitaev
        self.basis = "sto3g"

    def run(self):

        print("=======================")
        print("=== Starting run... ===")
        print("=======================")

        print("Start time: " + time.strftime("%b %d, %Y - %X"))
        start_time = time.time()

        print("Loading xyz file...")
        self.g = self.load_xyz(self.geometry_filename)
        self.driver = PySCFDriver(atom=self.convert_g_to_str(self.g),
                         unit=UnitsType.ANGSTROM,
                         charge=self.charge,
                         spin=self.spin,
                         max_cycle=5000,
                         max_memory=1024 * 128,
                         basis=self.basis
                         )

        def load_molecule(filename):
            if os.path.exists(filename):
                print(f"Found {filename}. Loading...")
                molecule = QMolecule(filename)
                molecule.load()
            else:
                # Regenerate
                print(f"Couldn't find {filename}. Regenerating...")
                molecule = self.driver.run()
                molecule.save(filename)
            return molecule

        print("Getting molecule...")
        self.molecule = load_molecule(self.filename)

        # ----- Generate circuit -----

        print("Constructing operator...")
        fermop, energy_shift = FermionicOperator.construct_operator(self.molecule, self.freeze_list, self.remove_list)

        # Generate qubit op
        print("Mapping fermion operator to qubit operator...")
        qubitop = fermop.mapping(self.map_type)

        # ----- Compile results -----

        end_time = time.time()

        results_dir = "results/" + time.strftime("%b-%d-%Y-%X")
        print("Saving results to " + results_dir)
        pathlib.Path(results_dir).mkdir(parents=True, exist_ok=True)

        # Write notes file
        with open(os.path.join(results_dir, "notes.txt"), "w") as f:
            f.write("==== RUN NOTES ====\n")
            f.write("Filename: " + str(self.filename) + "\n")
            f.write("Charge: " + str(self.charge) + "\n")
            f.write("Spin: " + str(self.spin) + "\n")
            f.write("Freeze list: " + str(self.freeze_list) + "\n")
            f.write("Remove list: " + str(self.remove_list) + "\n")
            f.write("Map type: " + str(self.map_type) + "\n")
            f.write("G: " + str(self.g) + "\n")
            f.write("Basis: " + str(self.basis) + "\n")
            f.write("\n")
            f.write("==== RUN DETAILS ====\n")
            f.write("Started: " + time.ctime(start_time) + "\n")
            f.write("Ended: " + time.ctime(end_time) + "\n")
            f.write("Run duration (in seconds): " + str(end_time - start_time) + "\n")
        
        # Save operator
        qubitop.to_file(os.path.join(results_dir, "qubitop.bin"))
        with open(os.path.join(results_dir, "qubitop.txt"), "w") as f:
            f.write(qubitop.print_details())

    # ---- Utility functions ----

    def convert_g_to_str(self, g):
        results = []
        for atom in g:
            [element, xyz] = atom
            [x, y, z] = xyz
            results.append(f"{element} {x} {y} {z}")
    
        return "; ".join(results)

    
    def load_xyz(self, xyz_fname):
        """
        Load a .xyz file into PySCF mol format
    
        """
    
        with open(xyz_fname, "r") as f:
            lines = f.readlines()
    
        # Remove header
        lines = lines[2:]
        output = []
        for line in lines:
            elems = line.split()
    
            atom = elems[0].capitalize()
            x = float(elems[1])
            y = float(elems[2])
            z = float(elems[3])
    
            output.append([atom, [x, y, z]])
    
        return output

if __name__ == "__main__":
    GenQubitOpLih().run()
