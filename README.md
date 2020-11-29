# QOSF-FeMoco2020

![FeMoco residue from 3U7Q (labeled ICS)](images/ics.png "FeMoco residue from 3U7Q (labeled ICS)")

Project for the QOSF mentorship program investigating the FeMo cofactor of nitrogenase

# Work Log

* (11/15/2020) Extracted ICS from [3U7Q](https://www.rcsb.org/3d-view/3U7Q/1) (see molecules/ICS.pdb)
    * This protein was found from [this paper](https://pubs.rsc.org/en/content/articlelanding/2019/cp/c8cp06930a#!divAbstract) and was chosen for being from a high-resolution x-ray diffraction (1.0 A) experiment and contained the correct atom (carbon) in the central cage.
* (11/15/2020) Added ICS visualization from JMOL (see images/ics.png)
* (11/15/2020) Added initial script to generate the Hamiltonian (see `gen_hamiltonian.py`)
* (11/19/2020) Added results from initial attempt at evaluating the Hamiltonian (SCF Convergence Issue)
* (11/25/2020) Added `casscf_test.py` to test CASSCF with ROHF. Added converged checkpoint file (`rohf.chk`).
* (11/25/2020) Ran STO-3G, S=3/2, Q=-1 ROHF (results in npy branch)

# General notes

## Generating the Hamiltonian

* Current Hamiltonian needs ~230-ish qubits
* Too big for most NISQ computers, but can run on D-Wave if the Hamiltonian is quadratic

## Thorneley-Lowe Model

Molybdenum nitrogenase performs the following reaction:

```
N2 + 8 H+ + 8 e− + 16 MgATP → 2 NH3 + H2 + 16 MgADP + 16 Pi
```

The full mechanism for this is unknown, but an experimentally-supported schematic is given by the [Lowe-Thorneley kinetic model](https://en.wikipedia.org/wiki/Nitrogenase#Lowe-Thorneley_kinetic_model):

![Lowe-Thorneley kinetic model](https://upload.wikimedia.org/wikipedia/en/a/a0/Lowe-Thorneley_Kinetic_Model.jpg)

E0 is the resting state for FeMoco.  E4 is the Janus state which is ready to accept N2 to produce ammonia.  A lot of focus is on the structure of the E4 state.

## Papers

TODO
