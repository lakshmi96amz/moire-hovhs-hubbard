#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt


def sym_gap_sweep(U_values):
    """
    For each U, return the symmetric-state gap (E_sym - E0) and the
    ground-space valley polarization magnitude |T_z|.
    """
    sym_gaps = np.zeros(len(U_values))
    abs_pol  = np.zeros(len(U_values))
    Pmat = valley_pol_matrix(sector_states, N_b)
    for i, U in enumerate(U_values):
        qop_U = to_qubit_op(H0_band, n, params, U)
        H_sec = sector_hamiltonian(qop_U, sector_states)
        E0, deg, Peigs = ground_space_polarization(H_sec, Pmat)
        E_sym = symmetric_state_energy(H_sec, sector_states, index_of,
                                       f1, f2, mK, mKp)
        sym_gaps[i] = E_sym - E0
        abs_pol[i]  = np.max(np.abs(Peigs))
    return sym_gaps, abs_pol


def plot_sym_gap(U_values, sym_gaps, abs_pol=None, savepath=None):
    """Plot the order parameter sym_gap vs U (and optionally |T_z|)."""
    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax1.plot(U_values, sym_gaps, "o-", color="C0",
             label=r"$E_{\mathrm{sym}} - E_0$")
    ax1.set_xlabel(r"$U$  (units of $v_F/a$)")
    ax1.set_ylabel(r"symmetric-state gap  $E_{\mathrm{sym}}-E_0$", color="C0")
    ax1.tick_params(axis="y", labelcolor="C0")

    if abs_pol is not None:
        ax2 = ax1.twinx()
        ax2.plot(U_values, abs_pol, "s--", color="C3", alpha=0.7,
                 label=r"$|T_z|$")
        ax2.set_ylabel(r"ground-space $|T_z|$", color="C3")
        ax2.tick_params(axis="y", labelcolor="C3")
        ax2.set_ylim(-0.05, 1.1)

    ax1.set_title("Valley-polarization order parameter vs U")
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=200, bbox_inches="tight")
    plt.show()
    return fig


if __name__ == "__main__":
    # fine sweep from 0 upward; keep the top below the manifold-to-band gap
    U_values = np.linspace(0.0, 1.5, 61)
    sym_gaps, abs_pol = sym_gap_sweep(U_values)

    print(f"{'U':>8} {'sym_gap':>12} {'|Tz|':>8}")
    for U, g, p in zip(U_values, sym_gaps, abs_pol):
        print(f"{U:8.4f} {g:12.6f} {p:8.4f}")

    # The gap is > 0 already at the smallest nonzero U: weak-coupling
    # instability characteristic of the HOVHS. No finite U_c is defined.
    print("\nNote: sym_gap lifts off zero immediately -> order at arbitrarily "
          "weak U (HOVHS-driven). No tolerance-independent U_c exists.")

    plot_sym_gap(U_values, sym_gaps, abs_pol, savepath="sym_gap_vs_U.png")

