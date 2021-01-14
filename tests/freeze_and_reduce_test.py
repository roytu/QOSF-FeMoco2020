
import sys
import os
import pathlib

import numpy as np

import logging
logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)
#logging.getLogger().setLevel(logging.INFO)


class TestFreezeAndReduce(object):
    def init(self):
        self.filename = "hdf5_files/lih_sto3g_(0,0).hdf5"
        self.freeze_list = range(0, 1)
        self.remove_list = range(5, 6)
        #self.freeze_list = range(184)
        #self.remove_list = range(190, 239)
        self.map_type = "parity" # parity, jordan_wigner, or bravyi_kitaev
        self.g = self.load_xyz("molecules/LiH.xyz")
 
    def test_real_qiskit(self):

        # ----- Import libraries and initialize -----

        from qiskit.chemistry.drivers import PySCFDriver, UnitsType
        from qiskit.chemistry import QMolecule, FermionicOperator
        from qiskit.chemistry.core import Hamiltonian, TransformationType, QubitMappingType
        from qiskit.chemistry.transformations import FermionicTransformation

        self.init()

        self.driver = PySCFDriver(atom=self.convert_g_to_str(self.g),
                         unit=UnitsType.ANGSTROM,
                         charge=0,
                         spin=0,
                         max_cycle=5000,
                         max_memory=1024 * 128,
                         basis='sto3g'
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

        self.molecule = load_molecule(self.filename)

        # ----- Perform test -----

        # Build qubit operator (full)
        fermop = FermionicOperator(
                h1=self.molecule.one_body_integrals,
                h2=self.molecule.two_body_integrals
        )
        
        # Freeze
        fermop, energy_shift = fermop.fermion_mode_freezing(self.freeze_list)
        
        # Remove
        fermop = fermop.fermion_mode_elimination(self.remove_list)
        
        # Generate qubit op
        qubitop = fermop.mapping('parity')
        #qubitop = Z2Symmetries.two_qubit_reduction(qubitop, num_particles)

        # Generate results
        qubit_op_str = qubitop.print_details()
        
        # Write results to file
        results_dir = "results/old"
        output_file = "qubitop1.txt"
        
        pathlib.Path(results_dir).mkdir(parents=True, exist_ok=True)
        with open(os.path.join(results_dir, output_file), "w") as f:
            f.write(qubit_op_str)
    
    def test_modified_qiskit(self):
        # ----- Import libraries and initialize -----

        sys.path.insert(0, os.getcwd())
        sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), "lib")))

        from lib.qiskit.chemistry.drivers import PySCFDriver, UnitsType
        from lib.qiskit.chemistry import QMolecule, FermionicOperator
        from lib.qiskit.chemistry.core import Hamiltonian, TransformationType, QubitMappingType
        from lib.qiskit.chemistry.transformations import FermionicTransformation

        # ----- Perform test -----

        self.init()

        self.driver = PySCFDriver(atom=self.convert_g_to_str(self.g),
                         unit=UnitsType.ANGSTROM,
                         charge=0,
                         spin=0,
                         max_cycle=5000,
                         max_memory=1024 * 128,
                         basis='sto3g'
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

        self.molecule = load_molecule(self.filename)

        # ----- Perform test -----

        fermop, energy_shift = FermionicOperator.construct_operator(self.molecule, self.freeze_list, self.remove_list)
        
        # Generate qubit op
        qubitop = fermop.mapping('parity')
        
        # Save operator
        #qubitop.to_file("results/new/femoco.qubitop")

        # Generate results
        qubit_op_str = qubitop.print_details()

        # Write results to file
        results_dir = "results/new"
        output_file = "qubitop1.txt"
        
        pathlib.Path(results_dir).mkdir(parents=True, exist_ok=True)
        with open(os.path.join(results_dir, output_file), "w") as f:
            f.write(qubit_op_str)

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

