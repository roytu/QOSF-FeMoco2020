
from qiskit.chemistry import FermionicOperator
from qiskit.chemistry.drivers import PySCFDriver, UnitsType
from qiskit.aqua.operators import Z2Symmetries
from util.mol import load_xyz

import logging
logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

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

# --- Create PySCF Driver and run ---

driver = PySCFDriver(atom=convert_g_to_str(g),
                     unit=UnitsType.ANGSTROM,
                     charge=-1,
                     spin=3,
                     max_cycle=5000,
                     max_memory=1024 * 128,
                     basis='sto3g'
                     )
                     
#driver = PySCFDriver(atom='H .0 .0 .0; H .0 .0 0.735',
#                     unit=UnitsType.ANGSTROM,
#                     basis='sto3g')
molecule = driver.run()
molecule.save("hdf5_files/femoco_sto3g_(-1,3).hdf5")

import pdb; pdb.set_trace()

# --- REST OF CODE IS RIPPED FROM SOMEWHERE ELSE
# --- so it is irrelevant. Feel free to edit this

num_particles = molecule.num_alpha + molecule.num_beta
num_spin_orbitals = molecule.num_orbitals * 2

# Build the qubit operator, which is the input to the VQE algorithm in Aqua
ferm_op = FermionicOperator(h1=molecule.one_body_integrals, h2=molecule.two_body_integrals)
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

