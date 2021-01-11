
import os
import sys
import numpy as np
from qiskit.chemistry import FermionicOperator
from qiskit.chemistry import QMolecule
from qiskit.chemistry.drivers import PySCFDriver, UnitsType
from qiskit.aqua.operators import Z2Symmetries
from util.mol import load_xyz
from qiskit.chemistry.core import Hamiltonian, TransformationType, QubitMappingType
from qiskit.chemistry.transformations import FermionicTransformation

import logging
logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)
#logging.getLogger().setLevel(logging.INFO)

# Use PySCF, a classical computational chemistry software
# package, to compute the one-body and two-body integrals in
# molecular-orbital basis, necessary to form the Fermionic operator

# --- Get molecule in string format (e.g. H .0 .0 1.0) ---
def convert_g_to_str(g):
    results = []
    for atom in g:
        [element, xyz] = atom
        [x, y, z] = xyz
        results.append(f"{element} {x} {y} {z}")

    return "; ".join(results)

g = load_xyz("molecules/ICS.xyz")
#g = load_xyz("molecules/ICS-min.xyz")
#g = load_xyz("molecules/LiH.xyz")

# --- Create PySCF Driver and run ---

driver = PySCFDriver(atom=convert_g_to_str(g),
                     unit=UnitsType.ANGSTROM,
                     charge=-1,
                     spin=3,
                     #charge=0,
                     #spin=0,
                     max_cycle=5000,
                     max_memory=1024 * 128,
                     basis='sto3g'
                     )
#driver = PySCFDriver(atom='H .0 .0 .0; H .0 .0 0.735',
#                     unit=UnitsType.ANGSTROM,
#                     basis='sto3g')
                     
FILENAME = "hdf5_files/femoco_sto3g_(-1,3).hdf5"
#FILENAME = "hdf5_files/femoco_nosulfer_sto3g_(-1,3).hdf5"
#FILENAME = "hdf5_files/lih_sto3g_(0,0).hdf5"
#FILENAME = "hdf5_files/hydrogen.hdf5"

if os.path.exists(FILENAME):
    print(f"Found {FILENAME}. Loading...")
    molecule = QMolecule(FILENAME)
    molecule.load()
else:
    # Regenerate
    print(f"Couldn't find {FILENAME}. Regenerating...")
    molecule = driver.run()
    molecule.save(FILENAME)

#print("Loading 1 body integrals...")
#one_body_integrals = molecule.one_body_integrals
##np.save("results/old/lih_1d.npy", one_body_integrals)
#np.save("results/new/lih_1d.npy", one_body_integrals)
#
#print("Loading 2 body integrals...")
#two_body_integrals = molecule.two_body_integrals
#np.save("results/new/lih_2d.npy", two_body_integrals)

#print("Loading FermionicOperator...")
#ferm_op = FermionicOperator(
#        h1=one_body_integrals,
#        h2=two_body_integrals
#)
#print("Loaded FermionicOperator")

#num_particles = molecule.num_alpha + molecule.num_beta
#num_spin_orbitals = molecule.num_orbitals * 2

#core = FermionicTransformation(
#    transformation=TransformationType.FULL,
#    qubit_mapping=QubitMappingType.PARITY, 
#    two_qubit_reduction=False,
#    orbital_reduction=range(0, 1),
#    freeze_core=True)
#qubit_op, aux_ops = core.transform(driver)
#sys.exit(0)

print("Start Hamiltonian Reduction Procedure")
#freeze_list = range(0, 1)
#remove_list = range(5, 6)
freeze_list = range(184)
remove_list = range(190, 239)
map_type = 'parity' # parity, jordan_wigner, or bravyi_kitaev

# Evaluate particle numbers
num_particles = molecule.num_alpha + molecule.num_beta
num_spin_orbitals = molecule.num_orbitals * 2 # TODO: treat UHF/ROHF later

if False:  # Old
	# Build qubit operator (full)
	fermop = FermionicOperator(h1 = molecule.one_body_integrals, h2 = molecule.two_body_integrals)

	# Freeze
	fermop, energy_shift = fermop.fermion_mode_freezing(freeze_list)

	# Remove
	fermop = fermop.fermion_mode_elimination(remove_list)

	# Update variables
	num_spin_orbitals -= len(freeze_list)
	num_particles -= len(freeze_list)
	num_spin_orbitals -= len(remove_list)
	
	# Generate qubit op
	qubitop = fermop.mapping('parity')
	qubitop = Z2Symmetries.two_qubit_reduction(qubitop, num_particles)

	qubit_op_str = qubitop.print_details()
	open("results/old/qubitop1.txt", "w").write(qubit_op_str)
else:  # New
	fermop, energy_shift = FermionicOperator.construct_operator(molecule, freeze_list, remove_list)

	# Update variables
	num_spin_orbitals -= len(freeze_list)
	num_particles -= len(freeze_list)
	num_spin_orbitals -= len(remove_list)

	# Generate qubit op
	import pdb; pdb.set_trace()
	qubitop = fermop.mapping('parity')
	qubitop = Z2Symmetries.two_qubit_reduction(qubitop, num_particles)
	import pdb; pdb.set_trace()

	# Save operator
	qubitop.to_file("results/new/femoco.qubitop")

	qubit_op_str = qubitop.print_details()
	#open("results/new/qubitop1.txt", "w").write(qubit_op_str)
	open("results/new/qubitop_femoco.txt", "w").write(qubit_op_str)

# ETC ETC OH GOD
sys.exit(0)
#qubitop.chop(10**-10)
## Calculate shift
#shift = energy_shift + nuclear_repulsion_energy


 
# Build the qubit operator, which is the input to the VQE algorithm in Aqua

map_type = 'PARITY'
qubit_op = ferm_op.mapping(map_type)
qubit_op = Z2Symmetries.two_qubit_reduction(qubit_op, num_particles)
num_qubits = qubit_op.num_qubits
sys.exit(0)



import pdb; pdb.set_trace()

core = Hamiltonian(
    transformation=TransformationType.FULL,
    qubit_mapping=QubitMappingType.PARITY, 
    two_qubit_reduction=False,
    orbital_reduction=range(0, 183),
    freeze_core=True)
qubit_op, aux_ops = core.run(molecule)

#import pdb; pdb.set_trace()


1/0

# HAMILTONIAN REDUCTION
# Specify which orbitals to freeze and remove.
# For example, in ICS.xyz,
# Orbitals with occupancy=2 and Mulliken orbital population near
# double occupancy can be readily frozen.
# Orbitals with occupancy=0 and Mulliken orbital population near
# zero occupancy can be readily removed (virtual orbitals).
# Following logic is in part adapted from the code in
# https://gist.github.com/0x6f736f646f/698ee32dde649ace70ad1152d276f748
# @PARAMS Specify orbital numbers and mapping type below.
print("Start Hamiltonian Reduction Procedure")
freeze_list = range(184)
remove_list = range(190, 239)
map_type = 'parity' # parity, jordan_wigner, or bravyi_kitaev
# Evaluate particle numbers
num_particles = molecule.num_alpha + molecule.num_beta
num_spin_orbitals = molecule.num_orbitals * 2 # TODO: treat UHF/ROHF later
# Build qubit operator (full)
fermop = FermionicOperator(h1 = molecule.one_body_integrals, h2 = molecule.two_body_integrals)
# Reduce
fermop, energy_shift = fermop.fermion_mode_freezing(freeze_list)
num_spin_orbitals -= len(freeze_list)
num_particles -= len(freeze_list)
fermop = fermop.fermion_mode_elimination(remove_list)
num_spin_orbitals -= len(remove_list)
qubitop = fermop.mapping('parity')
qubitOp = Z2Symmetries.two_qubit_reduction(qubitOp, num_particles)
qubitOp.chop(10**-10)
# Calculate shift
shift = energy_shift + nuclear_repulsion_energy


import pdb; pdb.set_trace()




# --- REST OF CODE IS RIPPED FROM SOMEWHERE ELSE
# --- so it is irrelevant. Feel free to edit this

num_particles = molecule.num_alpha + molecule.num_beta
num_spin_orbitals = molecule.num_orbitals * 2

# Build the qubit operator, which is the input to the VQE algorithm in Aqua
ferm_op = FermionicOperator(
        h1=molecule.one_body_integrals,
        h2=molecule.two_body_integrals
)
map_type = 'PARITY'
qubit_op = ferm_op.mapping(map_type)
qubit_op = Z2Symmetries.two_qubit_reduction(qubit_op, num_particles)
num_qubits = qubit_op.num_qubits

# setup a classical optimizer for VQE
from qiskit.aqua.components.optimizers import L_BFGS_B
optimizer = L_BFGS_B()

# setup the initial state for the variational form
from qiskit.chemistry.components.initial_states import HartreeFock
init_state = HartreeFock(num_spin_orbitals, num_particles)

# setup the variational form for VQE
from qiskit.circuit.library import TwoLocal
var_form = TwoLocal(num_qubits, ['ry', 'rz'], 'cz', initial_state=init_state)

# setup and run VQE
from qiskit.aqua.algorithms import VQE
algorithm = VQE(qubit_op, var_form, optimizer)

# set the backend for the quantum computation
from qiskit import Aer
backend = Aer.get_backend('statevector_simulator')

result = algorithm.run(backend)
print(result.eigenvalue.real)

