
import pyscf
from pyscf import gto, scf, mcscf

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

def build_femoco():
    """ Returns the PySCF Mol object for Femoco """

    g = load_xyz("../molecules/ICS.xyz")

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

