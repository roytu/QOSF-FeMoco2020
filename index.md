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

The purpose of this project was to calculate the ground state of various
derivatives of FeMoco, which is this molecule:

[TODO](image)

So what is this and why do we care? This is a cluster that appears in
nitrogenase, an enzyme responsible for converting gaseous nitrogen (N2) to
ammonia (NH3) in soil.  The ammonia is then used by plants to synthesize
chlorophyll. This process (called nitrogen fixation) is a major bottleneck for plant growth, and therefore
food production, and there are industrial processes that attempt to mimic this.
The Haber-Bosch process subjects gaseous nitrogen to high temperatures and
pressures to break the triple-bond, and produces the majority of the nitrogen
supply available to plants today.  However, generating this high-pressure
environment is energy-expensive, and the Haber process consumes about 1-2% of the
world's energy production.  What's /not/ understood is how nitrogenase can break
the triple-bond of N2 at atmospheric temperatures and pressures.

TODO add formula for Haber bosch
TODO add formula for Femoco catalyzed nitrogen fixation

This is the question: How does FeMoco *actually* catalyze nitrogen fixation, and
can we scale this process to replace Haber-Bosch?

## Background Info

To make this concrete -- we know that the reaction starts with N2 and FeMoco,
and somehow ends up producing NH3 and FeMoco:

TODO N2 + Femoco -> intermediates -> NH3 + Femoco image in Azobacter vinelandii

TODO We want to find what are the chemical intermediates
TODO Eyring rate equation
TODO Energy spectrum
TODO Lowest energy not necessarily ideal (reference Ian Dance paper)

## Project Planning

TODO

1. Run SCF on FeMoco
2. Generate Hamiltonian
3. Run VQE on hamiltonian

## Extracting the Molecule / Defining the Fock Space

The first step in analyzing the molecule is, naturally, to find the geometry for
it. The first structural models appeared in 1978[^cramer], with the six-atom
iron cage discovered in 1992. In 2002, a central atom within the cage was
discovered, which was widely believed to be nitrogen. It wasn't until 2011 when
we'd have the correct stoichiometry, when the central atom was determined
to be carbon [^spatzal].

All of this is to say, a lot of the literature and geometries for FeMoco are
incorrect, and it's important to find a structure from after 2011. Our analysis
uses [3U7Q](https://www.rcsb.org/structure/3u7q), a protein sample from
nitrogenase in Azotobacter vinelandii. The 3D view confirms that this has the
correct FeMoco cluster, labeled ICS 6496.

TODO add image of ICS 6496



## Generating the Hamiltonian

TODO

## Details on Qiskit / PySCF libraries

TODO

## Further Work

TODO

## Useful Papers

TODO


## Footnotes

[^qosf]: <https://qosf.org/qc_mentorship/>
[^vesselin]: <TODO>
[^spatzal]: <https://science.sciencemag.org/content/334/6058/940>
[^cramer]: . Cramer SP, Hodgson KO, Gillum WO, Mortenson LE. J Am Chem Soc.
[^bjornsson]: 1978;100:3398–3407. Bjornsson R, Neese F, Schrock RR, Einsle O, DeBeer S. The discovery of Mo(III) in FeMoco: reuniting enzyme and model chemistry. J Biol Inorg Chem. 2015;20(2):447-460. doi:10.1007/s00775-014-1230-6
