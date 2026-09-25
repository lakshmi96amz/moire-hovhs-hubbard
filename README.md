# moire-hovhs-hubbard
Exact-diagonalization and VQE study of an interaction-driven valley-polarized instability at a higher-order Van Hove singularity on a topological-insulator moir\'e surface state.

# Valley-polarized instability at a higher-order Van Hove singularity

This repository studies what a Hubbard interaction does to the higher-order Van Hove singularity (HOVHS) on the moir\'e surface state of a 3D topological insulator (Bi₂Te₃ parameters). Starting from the single-particle Chern-band model of [arXiv:2402.16772](https://arxiv.org/abs/2402.16772), which hosts a pair of cubic HOVHS at the $\pm K$ moir\'e valleys, we add an on-site Hubbard repulsion, project onto a two-valley (K, K′) patch, and solve the interacting problem by exact diagonalization in a fixed particle-number sector, validated with a variational quantum eigensolver (VQE).

## Result

At half-filling of the valley-degenerate HOVHS manifold, the Hubbard interaction drives **valley-Ising order**: the ground state spontaneously polarizes into a single valley (valley polarization $T_z = \pm 1$), and the two polarized states form a symmetry-protected degenerate doublet. The order is Ising rather than XY because the Hamiltonian conserves valley polarization (a U(1) valley symmetry) but not intervalley coherence. The instability sets in at arbitrarily weak coupling, the order parameter lifts off zero the moment U is turned on, as expected from the divergent density of states, is robust across multiple ($h_z$, V) points on the HOVHS line, and is independent of the sublattice polarization of the single-particle HOVHS states.


## Requirements

See `requirements.txt`. Core dependencies: `numpy`, `matplotlib`, `qiskit`, `qiskit-nature`, `qiskit-algorithms`. The Qiskit APIs used here (the `FermionicOp` constructor, the Jordan-Wigner mapper, and the VQE estimator) shift between releases, so the pinned versions in `requirements.txt` are the ones the code was run with.

## Reference

Built on L. Pullasseri and L. H. Santos, *Chern Bands with Higher-Order Van Hove Singularities on Topological Moiré Surface States*, [arXiv:2402.16772](https://arxiv.org/abs/2402.16772). This repository is the interacting follow-up to that single-particle work.
