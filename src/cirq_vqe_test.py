import math
import pickle
import itertools
import numpy as np
import openfermion
import openfermioncirq
import cirq
from time import time
from pyscf import gto, scf, fci
from copy import deepcopy
from openfermion import QubitOperator
from openfermion.hamiltonians import MolecularData
#from openfermionpyscf import run_pyscf
from run_pyscf_custom import run_pyscf
from openfermion.transforms import get_fermion_operator, bravyi_kitaev, jordan_wigner
from openfermion.utils import taper_off_qubits, commutator
from util.mol import load_xyz

def get_qubit_hamiltonian(g, basis, charge=0, spin=1, qubit_transf='jw'):

    ## Create OpenFermion molecule
    #mol = gto.Mole()
    #mol.atom = g
    #mol.basis = basis
    #mol.spin = spin
    #mol.charge = charge
    #mol.symmetry = False
    ##mol.max_memory = 1024
    #mol.build()

    multiplicity = spin + 1  # spin here is 2S ?
    mol = MolecularData(g, basis, multiplicity, charge)
    #mol.load()

    # Convert to PySCF molecule and run SCF
    print("Running run_pyscf...")
    print(f"Time: {time()}")
    print("=" * 20)
    mol = run_pyscf(mol)

    # Freeze some orbitals?
    occupied_indices = None
    active_indices = None

    #import pdb; pdb.set_trace()

    # Try:
    occupied_indices = range(183)
    active_indices = range(15)
    # when it stops lagging

    # Get Hamiltonian
    print("Running get_molecular_hamiltonian...")
    print(f"Time: {time()}")
    print("=" * 20)
    ham = mol.get_molecular_hamiltonian(
            occupied_indices=occupied_indices,
            active_indices=active_indices)

    print("Running get_fermion_operator...")
    print(f"Time: {time()}")
    print("=" * 20)
    hamf = get_fermion_operator(ham)

    print(f"Running {qubit_transf}...")
    print(f"Time: {time()}")
    print("=" * 20)

    if qubit_transf == 'bk':
        hamq = bravyi_kitaev(hamf)
    elif qubit_transf == 'jw':
        hamq = jordan_wigner(hamf)
    else:
        raise(ValueError(qubit_transf, 'Unknown transformation specified'))

    # Adapted from 10.5281/zenodo.2880550
    hamd = openfermion.get_diagonal_coulomb_hamiltonian(hamf, ignore_incompatible_terms=True)
    nqubits = openfermion.count_qubits(hamd)
    ansatz = openfermioncirq.SwapNetworkTrotterAnsatz(hamd, iterations=3)
    print(ansatz.circuit.to_text_diagram(transpose=True))
    
    nelectrons = 8 # TODO: CHECK WHICH ORBITALS ARE SAMPLED!!
    
    prep_circuit = cirq.Circuit(openfermioncirq.prepare_gaussian_state(ansatz.qubits, openfermion.QuadraticHamiltonian(hamd.one_body))) # TODO: Verify this line
    objective = openfermioncirq.HamiltonianObjective(hamd)
    study = openfermioncirq.VariationalStudy(name='hamil', ansatz=ansatz, objective=objective, preparation_circuit=prep_circuit)

    print("The energy of the default initial guess is {}.".format(
    study.value_of(ansatz.default_initial_params()))
)
    print()

    alg = openfermioncirq.optimization.ScipyOptimizationAlgorithm(kwargs={'method': 'COBYLA'}, uses_bounds=False)
    opt_params = openfermioncirq.optimization.OptimizationParams(algorith=alg, initial_guess=ansat.default_initial_params())
    result = study.optimize(opt_param)

    print("Optimized energy: {}".format(result.optimal_value))
    print()
    

    #return remove_complex(hamq)

def remove_complex(H : QubitOperator, tiny=1e-8):
    '''
    Removing near-zero complex coefficient
    '''
    real_h = QubitOperator.zero()
    for term, val in H.terms.items():
        if np.imag(val) < tiny:
            val = np.real(val)
        real_h += QubitOperator(term=term, coefficient=val)
    return real_h

if __name__ == "__main__":
    g = load_xyz("molecules/ICS.xyz")
    print(g)
    get_qubit_hamiltonian(g, "sto-3g", charge=-1, spin=3)
    #H = get_qubit_hamiltonian(g, "sto-3g", charge=-1, spin=3)
    #print(H)
