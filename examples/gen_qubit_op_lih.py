
from gen_qubit_op import GenQubitOp

class GenQubitOpLih(GenQubitOp):
    def __init__(self):
        super().__init__()

        self.filename = "hdf5_files/lih_sto3g_(0,0).hdf5"
        self.geometry_filename = "molecules/LiH.xyz"
        self.charge = 0
        self.spin = 0
        self.freeze_list = [2, 4]
        self.remove_list = [1, 9]
        self.map_type = "jordan_wigner" # parity, jordan_wigner, or bravyi_kitaev
        self.basis = "sto3g"

if __name__ == "__main__":
    GenQubitOpLih().run()
