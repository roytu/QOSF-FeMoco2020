# QOSF-FeMoco2020

![FeMoco residue from 3U7Q (labeled ICS)](images/ics.png "FeMoco residue from 3U7Q (labeled ICS)")

# What is this?

This is the project that Minsik Cho and I worked on during a mentorship program
hosted by the Quantum Open Software Foundation (<https://qosf.org>), over the
winter of the year 2020.  Our project attempts to apply the Variational Quantum
Eigensolver -- an algorithm to estimate energy levels of molecules -- to the
FeMo[^femoco] cofactor of nitrogenase.  This cluster is responsible for
nitrogenase's remarkable ability to break the triple bond of the `N_2` molecule
at atmospheric temperatures and pressures, and the ability to calculate energies
for chemical intermediates of the cluster is critical to understanding the
viable reaction pathways.

This is ongoing work, and may be in various states of disarray... be warned.

Visit our [Project Webpage](https://roytu.github.io/QOSF-FeMoco2020) for more
details.

# File Structure

* `examples/`: Example scripts
* `src/`: Larger proof-of-concept scripts (less stable)
* `lib/`: Modified qiskit code
* `tests/`: Unit tests. Run with `pytest tests/`
* `hdf5_files/`: SCF results. Load with `QMolecule.load()`
* `molecules/`: Molecular geometries (pulled from the [Protein Data Bank](https://www.rcsb.org/))

# Setup and Usage Instructions

1. Install [conda](https://docs.conda.io/en/latest/).
2. Clone this repository and `cd` to it.
3. Create the environment using:
    ```
    conda env create -f environment.yml
    ```
4. Activate the conda environment:
    ```
    conda activate qosf
    ```
5. Run the test suite:
    ```
    pytest tests/
    ```
6. Run the scripts from the repository root (e.g.):
    ```
    python examples/gen_qubit_op_lih.py
    python src/gen_hamiltonian_qiskit.py
    ```

# Special Thanks

Special thanks to [Vesselin](https://www.linkedin.com/in/vgg-consulting/) for providing guidance throughout this
project, [Michał Stęchły](https://www.mustythoughts.com/) for organizing the mentorship program,
and the rest of the folks at the [Quantum Open Software Foundation](qosf.org)
for their work in enabling free and open-source research for quantum computing.

# Contact Info

Questions? Reach us at:

* Roy Tu (kroytu [at] berkeley [dot] edu)
* [Minsik Cho](http://linkedin.com/in/chominsik)

[^femoco]: <https://en.wikipedia.org/wiki/FeMoco>
