
import pyscf
from pyscf import gto, scf, mcscf
from time import time

from util.mol import build_femoco

def casscf():
    print("Building molecule...")
    print(f"Time: {time()}")
    print("=" * 20)
    mol = build_femoco()

    print("Running ROHF...")
    print(f"Time: {time()}")
    print("=" * 20)
    mf = scf.ROHF(mol)
    mf.chkfile = '../chkfiles/rohf.chk'
    mf.init_guess = 'chkfile'

    # Just use initial guess
    mf.max_cycle = 0
    mf.verbose = 4
    mf.kernel()
    mf.analyze()

    #print("Running CASSCF...")
    #print(f"Time: {time()}")
    #print("=" * 20)

    #n_orbitals = 6
    #n_electrons = 6

    #mc = mcscf.CASSCF(mf, n_orbitals, n_electrons)
    #mc.chkfile = f'casscf-{n_orbitals}-{n_electrons}.chk'
    ##mc.init_guess = 'chkfile'
    #mc.verbose = 4
    #energy = mc.kernel()
    #print(energy)

if __name__ == "__main__":
    casscf()
