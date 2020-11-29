
import pyscf
from pyscf import gto, scf, mcscf
from gen_hamiltonian import load_xyz
from time import time

def build_femoco():
    """ Returns the PySCF Mol object for Femoco """

    g = load_xyz("molecules/ICS.xyz")

    mol = gto.Mole()
    mol.atom = g
    mol.basis = "sto-3g"
    mol.spin = 3
    mol.charge = -1
    mol.symmetry = True
    mol.max_memory = 1024 * 4
    mol.verbose = 4
    mol.build()

    return mol

def casscf():
    print("Building molecule...")
    print(f"Time: {time()}")
    print("=" * 20)
    mol = build_femoco()

    print("Running ROHF...")
    print(f"Time: {time()}")
    print("=" * 20)
    mf = scf.ROHF(mol)
    mf.chkfile = 'rohf.chk'
    mf.init_guess = 'chkfile'
    mf.max_cycle = 1000
    mf.verbose = 4
    mf.kernel()

    print("Running CASSCF...")
    print(f"Time: {time()}")
    print("=" * 20)

    n_orbitals = 6
    n_electrons = 6

    mc = mcscf.CASSCF(mf, n_orbitals, n_electrons)
    mc.chkfile = f'casscf-{n_orbitals}-{n_electrons}.chk'
    #mc.init_guess = 'chkfile'
    mc.verbose = 4
    energy = mc.kernel()
    print(energy)

if __name__ == "__main__":
    casscf()
