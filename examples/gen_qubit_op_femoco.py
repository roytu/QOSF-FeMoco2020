
from gen_qubit_op import GenQubitOp

class GenQubitOpFemoco(GenQubitOp):
    def __init__(self):
        super().__init__()

        self.filename = "hdf5_files/femoco_sto3g_(-1,3).hdf5"
        self.geometry_filename = "molecules/ICS.xyz"
        self.charge = -1
        self.spin = 3
        self.freeze_list = range(184)
        self.remove_list = range(190, 239)
        self.map_type = "jordan_wigner" # parity, jordan_wigner, or bravyi_kitaev
        self.basis = "sto3g"

if __name__ == "__main__":
    GenQubitOpFemoco().run()
