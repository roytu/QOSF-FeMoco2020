## Introduction

One of the most exciting applications of quantum computers is the ability to
simulate chemicals, at the layer of abstraction that nature itself operates.
Our current classical algorithms for molecular simulation -- while extremely
impressive -- will always be lacking, because a classical computer cannot
efficiently handle the exponential explosion in state space inherent in real
quantum systems.  For example, suppose we want to calculate the probability for
an electron at point X to move to point Y (i.e. `P(X -> Y)`). Nature's computer
calculates this by evaluating all probability amplitudes `P(X -> Z -> Y)` for
all intermediate points Z in space. A classical computer has no choice but to
iterate over every choice of Z, so the runtime will scale by `O(N)` where N is
the number of points in space we are considering.  This gets even worse when you
consider the second-order contributions (`P(X -> Z1 -> Z2 -> Y)`).

With a quantum computer, we can prepare a state with the electron at position X,
evolve it for a single time-step, and perform a measurement.  The time-evolution
performs all of the intermediate steps under the hood, and the measurement
collapses all of these states to return a new position Y.  By performing this
experiment many times, we can begin to reconstruct the probability distribution
`P(X -> Y)`.

This is not how classical or quantum simulations work in practice, but serves as
an illustration for why these calculations are so hard for classical computers
to deal with.

TODO add the following based on length:
- Rise of NISQ systems
- Rise of VQE
    - Talk about QPE?

--- 

This post serves to aggregate what we learned from a three-month mentorship
hosted by the Quantum Open Software Foundation [^qosf], through interactions
with our mentor (Vesselin [^vesselin]) and the other mentees. We are forever
grateful for the time and energy they have donated, for which none of this would
have been possible.

During the mentorship, we investigated the viability of simulating
[FeMoco](TODO), a cluster of molecules in nitrogenase.  The cluster is
responsible for the nitrogenase enzyme's remarkable ability to catalyze nitrogen
fixation, so there is significant interest in accurately simulating this
molecule.  We've structured this post to target someone who has /some/
background in quantum computation and computational chemistry, but may not have
a full grasp of the process for performing this simulation, and what roadblocks
arise.  It is more of a pedagogical journey than a textbook -- you should feel
free to skip sections as needed.  So without further ado...

## Problem Overview

TODO

## Background Info

### Quantum Chemistry using Classical Computers
Simulating quantum mechanical systems, such as molecules have been an important quest in the theoretical chemistry community for some time and numerous algorithms of different cost and underlying theories, are actively researched. For transition metal complexes, due to its complexity, a common choice is to use DFT (density functional theory). Starting from the many-body setup with multiple electrons, protons, and neutrons interacting in the target system, DFT operates on the electron density as a variable, reducing down some complexity of the problem. Based on this background, DFT has some record of treating some moderately-sized molecules chemists are interested in with reasonable cost and accuracy in predicting key chemical properties.

Even with these simplifications, the entire FeMoco complex is quite difficult to simulate with DFT using basis sets that are detailed enough and clever choice of the functionals. Because our aim is to pedagogically explore the steps as if we had the resources to complete the simulation on conventional and quantum computers, simple *(and fictitious)* complex of Fe<sub>2</sub>S<sub>2</sub>, which comprise the core of the FeMoco was used as a model system. Using Orca[^orca] and Avogadro[^avogadro], starting from initial geometries deposited from crystallography, geometry optimization can be performed using the functionals and basis sets (TPSSh and def2-TZVPP in our run) and the frontier orbitals can be plotted. These simulations on classical computer provide us with a baseline of where things are before we explore the quantum computing side of this tale.

[//]: # "TODO Add diagrams/visuals from DFT run"

### The Variational Principle

[//]: # "TODO Add a few words about variational principle"

### Variational Quantum Eigensolver

[//]: # "TODO Add a few words about VQE"

## Project Planning

TODO

## Generating the Hamiltonian

TODO

## Details on Qiskit / PySCF libraries

TODO

## Further Work

TODO


## Footnotes

[^qosf]: <https://qosf.org/qc_mentorship/>
[^vesselin]: <TODO>
[^orca]: <https://orcaforum.kofo.mpg.de/>
[^avogadro]: <https://avogadro.cc/>
