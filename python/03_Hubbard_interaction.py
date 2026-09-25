#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.mappers   import JordanWignerMapper

def _mode(ik, r):
    """Flat mode index for patch point ik and band r."""
    return ik * N_b + r

def build_density_op(n, params, i_out, i_in, sigma):
    """
    Band-projected density operator matrix block connecting patch point
    i_in (annihilation) to patch point i_out (creation), for spin sigma:

        rho_sigma(i_out <- i_in)
           = sum_{a,b} M_sigma(k_out, a ; k_in, b) gamma^dag_{i_out,a} gamma_{i_in,b}

    Returned as a FermionicOp on N_modes spin-orbitals.
    """
    k_out = Kpatch[i_out]
    k_in  = Kpatch[i_in]
    data = {}
    for a, b in product(range(N_b), repeat=2):
        # M_sigma(k_out, a ; k_in, b) = sum_G conj(U^sigma_{G,a}(k_out)) U^sigma_{G,b}(k_in)
        coeff = get_overlap(n, params, a, b, k_out, k_in, sigma)
        if abs(coeff) < 1e-14:
            continue
        p = _mode(i_out, a)   # creation
        qd = _mode(i_in, b)   # annihilation
        key = f"+_{p} -_{qd}"
        data[key] = data.get(key, 0.0) + coeff
    return FermionicOp(data, num_spin_orbitals=N_modes)

def build_H_int(n, params, U):
    """
    Hubbard interaction as a FermionicOp, built from density-operator products.
    Hermitian by construction; includes normal-ordering one-body pieces.
    """
    H_int = FermionicOp({}, num_spin_orbitals=N_modes)

    for ik, ikp, ikq in product(range(N_k), repeat=3):
        k     = Kpatch[ik]
        kp    = Kpatch[ikp]
        kq    = Kpatch[ikq]                     # k + q
        kpmq  = [kp[0] - (kq[0] - k[0]),
                 kp[1] - (kq[1] - k[1])]        # k' - q  (raw)
        ikpq  = find_patch_index_folded(kpmq, Kpatch, b1, b2)
        if ikpq is None:
            continue

        # up-density block: creation at (k+q)=ikq, annihilation at k=ik
        rho_up   = build_density_op(n, params, ikq,  ik,  0)
        # down-density block: creation at (k'-q)=ikpq, annihilation at k'=ikp
        rho_down = build_density_op(n, params, ikpq, ikp, 1)

        # operator product; library normal-orders it
        prod = rho_up @ rho_down
        H_int = H_int + 0.5 * (prod + prod.adjoint())

    # prefactor
    H_int = (U / N_k) * H_int
    return H_int.simplify()

#Compuing the non-interacting part of the Hamiltonian in the band basis
def build_H0_band(n, N_b, params):
    dim = N_k * N_b;
    H0 = np.zeros((dim,dim));
    for ik in range(N_k):
        for r in range(N_b):
            m = ik*N_b + r
            Ham = build_hamiltonian(n,Kpatch[ik],params);
            H0[m,m] = H0[m,m] = np.linalg.eigvalsh(Ham)[bands[r]]

    return H0


def build_H0_fermionic(H0_band):
    """One-body diagonal H0 as a FermionicOp."""
    data = {}
    for m in range(N_modes):
        e = float(np.real(H0_band[m, m]))
        if abs(e) > 1e-14:
            data[f"+_{m} -_{m}"] = e
    return FermionicOp(data, num_spin_orbitals=N_modes)

def to_qubit_op(H0_band, n, params, U):
    """Full H = H0 + H_int, Jordan-Wigner mapped to a SparsePauliOp."""
    H0_op  = build_H0_fermionic(H0_band)
    Hint   = build_H_int(n, params, U)
    H_full = (H0_op + Hint).simplify()
    return JordanWignerMapper().map(H_full)

