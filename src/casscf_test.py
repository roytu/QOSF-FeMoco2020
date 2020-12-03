
import os
import pyscf
from pyscf import gto, scf, mcscf
from time import time

from util.mol import build_femoco

def casscf():
    print("Building molecule...")
    print(f"Time: {time()}")
    print("=" * 20)
    mol = build_femoco()

    ####### RUN ROHF #######
    mf = scf.ROHF(mol)

    # Try to load ROHF object if it exists
    chkfile_path = "../chkfiles/rohf.chk"
    if os.path.exists(chkfile_path):
        print(f"Found {chkfile_path}. Loaded.")
        print("Running ROHF...")
        print(f"Time: {time()}")
        print("=" * 20)
        mf.chkfile = chkfile_path
        mf.init_guess = "chkfile"
        mf.max_cycle = 0
    else:
        print(f"No {chkfile_path} found.")
        print("Running ROHF...")
        print(f"Time: {time()}")
        print("=" * 20)
        mf.chkfile = chkfile_path
        mf.init_guess = "chkfile"

    mf.scf()
    mf.analyze()

    ####### RUN CASSCF #######
    n_orbitals = 6
    n_electrons = 7

    mc = mcscf.CASSCF(mf, n_orbitals, n_electrons)

    chkfile_path = f'../chkfiles/casscf-{n_orbitals}-{n_electrons}.chk'

    # Try to load CASSCF object if it exists
    if os.path.exists(chkfile_path):
        print(f"Found {chkfile_path}. Loaded.")
        print("Running CASSCF...")
        print(f"Time: {time()}")
        print("=" * 20)
        mc.chkfile = chkfile_path
        mc.init_guess = "chkfile"
        mc.max_cycle = 0
    else:
        print(f"No {chkfile_path} found.")
        print("Running CASSCF...")
        print(f"Time: {time()}")
        print("=" * 20)
        mc.chkfile = chkfile_path
        mc.init_guess = "chkfile"

    mc.kernel()
    mc.analyze()

if __name__ == "__main__":
    casscf()
