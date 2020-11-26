import math
import pickle
import itertools
import numpy as np
import openfermion
from pyscf import gto, scf, fci
from copy import deepcopy
from openfermion import QubitOperator
from openfermion.hamiltonians import MolecularData
from openfermionpyscf import run_pyscf
from openfermion.transforms import get_fermion_operator, bravyi_kitaev, jordan_wigner
from openfermion.utils import taper_off_qubits, commutator


def load_xyz(xyz_fname):
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

def get_qubit_hamiltonian(g, basis, charge=0, spin=1, qubit_transf='jw'):

    mol = gto.Mole()
    mol.atom = g
    mol.basis = basis
    mol.spin = spin
    mol.charge = charge
    mol.symmetry = True
    #mol.max_memory = 1024
    mol.build()

    print(mol)

    ham = mol.get_molecular_hamiltonian()
    print(ham)
    hamf = get_fermion_operator(ham)

    if qubit_transf == 'bk':
        hamq = bravyi_kitaev(hamf)
    elif qubit_transf == 'jw':
        hamq = jordan_wigner(hamf)
    else:
        raise(ValueError(qubit_transf, 'Unknown transformation specified'))

    return remove_complex(hamq)

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
    H = get_qubit_hamiltonian(g, "sto-3g", charge=-1, spin=3)
    print(H)
