#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np


def ising_classification(H_sec, sector_states, index_of, N_b, mK, mKp):
    """
    Return the ground-space pseudospin eigenvalues and the full-sector
    commutators [H, T_a], which together classify the order.
    """
    Pmat = valley_pol_matrix(sector_states, N_b)
    Tx, Ty = build_coherence_operators(sector_states, index_of, mK, mKp)

    Tz_eigs, _ = project_into_ground_space(H_sec, Pmat)
    Tx_eigs, _ = project_into_ground_space(H_sec, Tx)
    Ty_eigs, _ = project_into_ground_space(H_sec, Ty)

    commutators = {
        "[H, T_z]": commutator_norm(H_sec, Pmat),
        "[H, T_x]": commutator_norm(H_sec, Tx),
        "[H, T_y]": commutator_norm(H_sec, Ty),
    }
    return {
        "Tz_ground": Tz_eigs,
        "Tx_ground": Tx_eigs,
        "Ty_ground": Ty_eigs,
        "commutators": commutators,
    }


if __name__ == "__main__":
    res = ising_classification(H_sec_U, sector_states, index_of, N_b, mK, mKp)

    print("ground-space pseudospin eigenvalues:")
    print("  T_z:", np.round(res["Tz_ground"], 6))
    print("  T_x:", np.round(res["Tx_ground"], 6))
    print("  T_y:", np.round(res["Ty_ground"], 6))
    print("\nfull-sector commutators (conserved quantity <=> ~0):")
    for k, v in res["commutators"].items():
        print(f"  {k}: {v:.2e}")

    tz_c = res["commutators"]["[H, T_z]"]
    tx_c = res["commutators"]["[H, T_x]"]
    ty_c = res["commutators"]["[H, T_y]"]
    if tz_c < 1e-8 and tx_c > 1e-6 and ty_c > 1e-6:
        print("\n=> T_z conserved, T_x/T_y not: VALLEY-ISING order "
              "(discrete K/K' polarization, U(1) valley symmetry).")
    else:
        print("\n=> commutation pattern differs; re-examine the symmetry.")

