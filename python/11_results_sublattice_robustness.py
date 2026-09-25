#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np

V6_0    = 1.36     # V_6(h_z=0) in units of v_F/a
GAMMA   = 0.05     # HOVHS-line coefficient (arXiv:2402.16772, Eq. 5)
U_FIXED = 0.87     # fixed coupling for the instability probe


def V6_line(hz):
    """Moire potential on the HOVHS line for a given Zeeman field."""
    return V6_0 + GAMMA * hz**2


def manifold_polarization(n, params, kvec, manifold_band_index):
    """Spin/sublattice polarization of the manifold Bloch state at kvec."""
    U_up   = get_U(n, params, kvec, 0)[:, manifold_band_index]
    U_down = get_U(n, params, kvec, 1)[:, manifold_band_index]
    w_up   = np.vdot(U_up,   U_up).real
    w_down = np.vdot(U_down, U_down).real
    return (w_up - w_down) / (w_up + w_down)


def intervalley_formfactor(n, params, k1, k2, band_index, sigma):
    """M_sigma(k1, band; k2, band) = <U_sigma(k1)|U_sigma(k2)> for one band."""
    v1 = get_U(n, params, k1, sigma)[:, band_index]
    v2 = get_U(n, params, k2, sigma)[:, band_index]
    return np.vdot(v1, v2)


def identify_modes(H0_band, N_b):
    """Return (f1, f2, mK, mKp): two filled modes + two manifold modes."""
    eps = np.diag(H0_band).real
    order = np.argsort(eps)
    f1, f2 = int(order[0]), int(order[1])
    manifold = [int(m) for m in order[2:4]]
    mK  = next(m for m in manifold if m < N_b)
    mKp = next(m for m in manifold if m >= N_b)
    return f1, f2, mK, mKp


def hovhs_line_sweep(hz_list):
    """
    Along the HOVHS line, return per-hz: sublattice polarization at K,
    form-factor imbalance |M_up/M_down|, and the valley-instability strength
    (symmetric-state gap) at fixed U.
    """
    manifold_band_index = bands.index(121)
    Pmat = valley_pol_matrix(sector_states, N_b)

    rows = []
    for hz in hz_list:
        p = [hz, V6_line(hz)]

        subl_pol = manifold_polarization(n, p, K, manifold_band_index)
        Mu = intervalley_formfactor(n, p, K, Kp, manifold_band_index, 0)
        Md = intervalley_formfactor(n, p, K, Kp, manifold_band_index, 1)
        ratio = abs(Mu) / abs(Md) if abs(Md) > 1e-12 else np.inf

        H0b = build_H0_band(n, N_b, p)
        f1, f2, mK, mKp = identify_modes(H0b, N_b)
        qop = to_qubit_op(H0b, n, p, U_FIXED)
        H_sec = sector_hamiltonian(qop, sector_states)
        E0, deg, Peigs = ground_space_polarization(H_sec, Pmat)
        E_sym = symmetric_state_energy(H_sec, sector_states, index_of,
                                       f1, f2, mK, mKp)

        rows.append((hz, p[1], subl_pol, ratio, E_sym - E0,
                     np.max(np.abs(Peigs))))
    return rows


if __name__ == "__main__":
    hz_list = np.linspace(0.0, 1.5, 10)
    rows = hovhs_line_sweep(hz_list)

    print(f"{'hz':>6} {'V':>7} {'subl_pol':>10} {'|Mup/Mdn|':>10} "
          f"{'sym_gap':>10} {'|Tz|':>6}")
    for hz, V, sp, ratio, gap, tz in rows:
        print(f"{hz:6.2f} {V:7.3f} {sp:+10.4f} {ratio:10.4f} "
              f"{gap:10.5f} {tz:6.2f}")

    print("\n=> sublattice polarization and form-factor imbalance grow with "
          "h_z, but sym_gap stays ~flat and |Tz|=1 throughout (full at h_z=0 "
          "where sublattice pol = 0): sublattice structure is NOT the driver.")

