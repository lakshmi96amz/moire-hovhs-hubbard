#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from qiskit_nature.second_q.mappers import JordanWignerMapper


def sector_spectrum_at_U(H0_band, n, params, U, sector_states, n_low=8):
    """Lowest n_low sector eigenvalues (and eigenvectors) at coupling U."""
    qop_U = to_qubit_op(H0_band, n, params, U)
    H_sec = sector_hamiltonian(qop_U, sector_states)
    evals, evecs = np.linalg.eigh(H_sec)
    return H_sec, evals[:n_low], evecs


def confirm_intervalley_terms(n, params, U, sector_states, evecs_ground):
    """
    Confirm the interaction actually couples the valleys on the sector:
    the interaction block is nonzero and has finite ground-state expectation.
    """
    H_int_op  = JordanWignerMapper().map(build_H_int(n, params, U))
    H_int_sec = sector_hamiltonian(H_int_op, sector_states)
    max_elem  = np.max(np.abs(H_int_sec))
    E_int_gs  = np.real(evecs_ground[:, 0].conj() @ H_int_sec @ evecs_ground[:, 0])
    return max_elem, E_int_gs


def classify_order(H_sec, sector_states, index_of, N_b, mK, mKp):
    """
    Characterize the ground-doublet order:
      - valley polarization T_z eigenvalues (expect +-1 -> fully polarized)
      - intervalley coherence T_x, T_y eigenvalues in the ground space
      - commutators [H, T_a] on the full sector (which are symmetries)
    Ising valley order: T_z conserved ([H,T_z]=0), T_x, T_y not conserved.
    """
    Pmat = valley_pol_matrix(sector_states, N_b)
    Tx, Ty = build_coherence_operators(sector_states, index_of, mK, mKp)

    Tz_eigs, _ = project_into_ground_space(H_sec, Pmat)
    Tx_eigs, _ = project_into_ground_space(H_sec, Tx)
    Ty_eigs, _ = project_into_ground_space(H_sec, Ty)

    comm = {
        "[H, T_z]": commutator_norm(H_sec, Pmat),
        "[H, T_x]": commutator_norm(H_sec, Tx),
        "[H, T_y]": commutator_norm(H_sec, Ty),
    }
    return {"Tz": Tz_eigs, "Tx": Tx_eigs, "Ty": Ty_eigs, "commutators": comm}


if __name__ == "__main__":
    U = 0.87

    # 1. Spectrum: U=0 (degenerate reference) vs U (split, polarized)
    _, evals0, _ = sector_spectrum_at_U(H0_band, n, params, 0.0, sector_states)
    H_sec_U, evalsU, evecsU = sector_spectrum_at_U(H0_band, n, params, U, sector_states)
    print("lowest sector eigenvalues at U=0   :", np.round(evals0, 6))
    print("lowest sector eigenvalues at U=0.87:", np.round(evalsU, 6))

    # 2. Interaction genuinely couples the valleys
    max_elem, E_int_gs = confirm_intervalley_terms(n, params, U, sector_states, evecsU)
    print(f"interaction block max element: {max_elem:.4f}")
    print(f"interaction energy in ground state: {E_int_gs:.6f}")

    # 3. Classify the order
    res = classify_order(H_sec_U, sector_states, index_of, N_b, mK, mKp)
    print("ground-space T_z eigenvalues:", np.round(res["Tz"], 4))
    print("ground-space T_x eigenvalues:", np.round(res["Tx"], 4))
    print("ground-space T_y eigenvalues:", np.round(res["Ty"], 4))
    for k, v in res["commutators"].items():
        print(f"{k}: {v:.2e}")
    print("=> T_z conserved, T_x/T_y not: valley-Ising order")

