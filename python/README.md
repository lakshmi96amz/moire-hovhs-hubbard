## Layout of the python files
00_Parameters.py: All physical and numerical parameters defined including the plane-wave cutoff, Hamiltonian parameters, the two-valley patch points and the chosen bands.

01_Single_particle_hamiltonian.py: Constructs $H = v_F(k\times\sigma)_z + V(r)\sigma_0 + h_z\sigma_z$ in the plane-wave basis and returns band energies. This is the non-interacting model of arXiv:2402.16772.

02_Band_basis.py: Projects the plane-wave eigenvectors onto the kept bands at each valley, computes the spin-projected overlap form factors $M_\sigma$, and handles momentum conservation with umklapp folding.

03_Hubbard_interaction.py: Builds the interaction as $H_{int} = (U/N_k)\sum_q \rho_\uparrow(q)\rho_\downarrow(-q)$ from the density operators, which is Hermitian by construction and includes the normal-ordering correction automatically. Assembles $H_0 + H_{int}$ and maps to qubits via Jordan-Wigner.

04_Validation.py: Validation tests for the band basis and the density operator construction.

05_Manifold_identification.py: Inspecting the single-particle spectrum to locate the HOVHS manifold and confirming the Fermi level sits at E_K.

06_FixedN_sector.py: Restricts to the half-filled manifold sector (N=3), where the Fermi level sits at the HOVHS energy $E_K$. Builds the sector basis and the sector-restricted Hamiltonian.

07_Order_parameters.py: Valley polarization $T_z$ and intervalley coherence $T_x, T_y$ as sector operators, projected into the ground space (basis-independent) to characterize the order.

08_results_instability.py: Turning on U shows the ground doublet becomes fully valley-polarized ($T_z = \pm 1$)

09_results_u_sweep.py: Sweeps the Hubbard U and tracks the valley-symmetric-state energy gap,$sym_gap(U) = E_sym(U) - E0(U)$, which is the order parameter for valley polarization.

10_results_ising_classification.py: H conserves $T_z$ but not $T_x, T_y$: the order is valley-Ising, a discrete K/K' population imbalance protected by U(1) valley symmetry.

11_results_sublattice_robustness.py: Along the HOVHS line, the valley instability is independent of the sublattice polarization of the manifold states.

12_results_vqe_validation.py: A particle-number-preserving VQE reproduces the exact-diagonalization ground energy and valley polarization in the N=3 sector.
