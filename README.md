# moire-hovhs-hubbard
Exact-diagonalization and VQE study of an interaction-driven valley-polarized instability at a higher-order Van Hove singularity on a topological-insulator moir\'e surface state.

# Valley-polarized instability at a higher-order Van Hove singularity

This repository studies what a Hubbard interaction does to the higher-order Van Hove singularity (HOVHS) on the moir\'e surface state of a 3D topological insulator (Bi₂Te₃ parameters). Starting from the single-particle Chern-band model of [arXiv:2402.16772](https://arxiv.org/abs/2402.16772), which hosts a pair of cubic HOVHS at the $\pm K$ moir\'e valleys, we add an on-site Hubbard repulsion, project onto a two-valley (K, K′) patch, and solve the interacting problem by exact diagonalization in a fixed particle-number sector, validated with a variational quantum eigensolver (VQE).

## Main result

At half-filling of the valley-degenerate HOVHS manifold, the Hubbard interaction drives **valley-Ising order**: the ground state spontaneously polarizes into a single valley (valley polarization $T_z = \pm 1$), and the two polarized states form a symmetry-protected degenerate doublet. The order is Ising rather than XY because the Hamiltonian conserves valley polarization (a U(1) valley symmetry) but not intervalley coherence. The instability sets in at arbitrarily weak coupling, the order parameter lifts off zero the moment U is turned on, as expected from the divergent density of states, is robust across multiple ($h_z$, V) points on the HOVHS line, and is independent of the sublattice polarization of the single-particle HOVHS states.

## Method

The interaction is built from band-projected density operators,

```
H_int = (U / N_k) Σ_q ρ_↑(q) ρ_↓(−q),
```

which is Hermitian by construction and generates the normal-ordering one-body correction automatically. Momentum conservation on the two-valley patch, including umklapp folding (using 3K ≡ 0 modulo the moiré reciprocal lattice), is handled explicitly. The many-body problem is solved by exact diagonalization in the fixed N = 3 sector, which places the Fermi level at the HOVHS energy $E_K$. A particle-number-penalized `EfficientSU2` VQE reproduces the exact ground-state energy, establishing that the instability is accessible to near-term quantum algorithms.

The order is characterized with three ground-space measurements: the valley polarization $T_z$ (diagonal), the intervalley coherence $T_x, T_y$ (off-diagonal), and the commutators $[H, T_a]$, which distinguish Ising from XY order. Projecting an operator into the two-fold degenerate ground space gives a basis-independent readout; the commutation with H, not the projected eigenvalues, is what fixes the universality class.

## Layout of the python files
00_Parameters.py: All physical and numerical parameters defined including the plane-wave cutoff, Hamiltonian parameters, the two-valley patch points and the chosen bands.
01_Single_particle_hamiltonian.py: Constructs $H = v_F(k\times\sigma)_z + V(r)\sigma_0 + h_z\sigma_z$ in the plane-wave basis and returns band energies. This is the non-interacting model of arXiv:2402.16772.
02_Band_basis.py: Projects the plane-wave eigenvectors onto the kept bands at each valley, computes the spin-projected overlap form factors $M_\sigma$, and handles momentum conservation with umklapp folding.
03_Hubbard_interaction.py: Builds the interaction as $H_{int} = (U/N_k)\sum_q \rho_\uparrow(q)\rho_\downarrow(-q)$ from the density operators, which is Hermitian by construction and includes the normal-ordering correction automatically. Assembles $H_0 + H_{int}$ and maps to qubits via Jordan-Wigner.
04_Validation.py: 
05_Manifold_identification.py: 
06_FixedN_sector.py: Restricts to the half-filled manifold sector (N=3), where the Fermi level sits at the HOVHS energy $E_K$. Builds the sector basis and the sector-restricted Hamiltonian.
07_Order_parameters.py: Valley polarization $T_z$ and intervalley coherence $T_x, T_y$ as sector operators, projected into the ground space (basis-independent) to characterize the order.
08_results_instability.py: Turning on U splits... shows the ground doublet becomes fully valley-polarized ($T_z = \pm 1$)
09_results_u_sweep.py: Sweeps the Hubbard U and tracks the valley-symmetric-state energy gap,$sym_gap(U) = E_sym(U) - E0(U)$, which is the order parameter for valley polarization.
10_results_ising_classification.py: H conserves $T_z$ but not $T_x, T_y$: the order is valley-Ising, a discrete K/K' population imbalance protected by U(1) valley symmetry.
11_results_sublattice_robustness.py: Along the HOVHS line, the valley instability is independent of the sublattice polarization of the manifold states.
12_results_vqe_validation.py: A particle-number-preserving VQE reproduces the exact-diagonalization ground energy and valley polarization in the N=3 sector.


## Requirements

See `requirements.txt`. Core dependencies: `numpy`, `matplotlib`, `qiskit`, `qiskit-nature`, `qiskit-algorithms`. The Qiskit APIs used here (the `FermionicOp` constructor, the Jordan-Wigner mapper, and the VQE estimator) shift between releases, so the pinned versions in `requirements.txt` are the ones the code was run with.

## Reference

Built on L. Pullasseri and L. H. Santos, *Chern Bands with Higher-Order Van Hove Singularities on Topological Moiré Surface States*, [arXiv:2402.16772](https://arxiv.org/abs/2402.16772). This repository is the interacting follow-up to that single-particle work.
