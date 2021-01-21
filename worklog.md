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
* (12/03/2020) Some ROHF runs and added Mulliken analysis results
* (12/03/2020) Added Fermionic Hamiltonian with 183 frozen orbitals, 30 active orbitals (`jw_ham`)
* (12/13/2020) Added Hamiltonian with 183 frozen orbitals, 15 active orbitals
* (12/21/2020) Added `hamiltonians/hamq_occ183_act5_jw_fin`, generated with:

```
hamf_occ183_act5 = get_fermion_operator(ham_occ183_act5)
hamq_occ183_act5_jw = jordan_wigner(hamf_occ183_act5)
hamq_occ183_act5_jw_fin = remove_complex(hamq_occ183_act5_jw)
x = str(hamq_occ183_act5_jw_fin)
x = hamq_occ183_act5_jw_fin
y = str(x)
open("hamq_occ183_act5_jw_fin", "w").write(y)
```
* (12/21/2020) Trying to rewrite code in qiskit/pyscf. Added `src/gen_hamiltonian_qiskit.py`

* (01/01/2020 - 01/11/2020) Rewriting qiskit (TODO back-update this)
* (01/11/2020) (develop-roy) Qiskit rewrite to reduce O((M+N)^4) freezing code to O(M^2 N), where

```
M = num. of non-frozen orbitals
N = num. of frozen orbitals
```

This should lead to a drastic reduction in memory and run-time requirements.  For FeMoco, we are trying values on the order of M=10, N=200.  This reduces the number of floating points from `(200+10)^4 = 1,944,810,000` to (10^2) * 200 = 20,000`

* (01/14/2020) (develop-roy) Added test code. Run `pytest tests/`

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

## Mulliken Analysis / Frozen Core Model

TODO

## QUBO

TODO
