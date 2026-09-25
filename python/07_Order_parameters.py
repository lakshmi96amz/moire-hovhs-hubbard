#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np


#Valley polarization T_z (diagonal on the sector)
def valley_pol_matrix(sector_states, N_b):
    """T_z = N_K - N_K' as a diagonal matrix on the sector basis."""
    diag = np.zeros(len(sector_states))
    for i, s in enumerate(sector_states):
        nK  = sum((s >> m) & 1 for m in range(0, N_b))
        nKp = sum((s >> m) & 1 for m in range(N_b, 2 * N_b))
        diag[i] = nK - nKp
    return np.diag(diag)


#Intervalley coherence T_x, T_y (off-diagonal hopping)
def jw_sign(state, m):
    """Jordan-Wigner sign for acting on mode m: (-1)^(occupied modes < m)."""
    mask = (1 << m) - 1
    return -1.0 if bin(state & mask).count("1") % 2 else 1.0


def apply_cdag_c(state, p, q):
    """
    Apply c^dag_p c_q to an occupation-string `state`.
    Returns (new_state, amplitude), or (None, 0.0) if annihilated.
    Mode m occupied iff bit m set; Jordan-Wigner signs included.
    """
    if not ((state >> q) & 1):            # c_q needs q occupied
        return None, 0.0
    sgn = jw_sign(state, q)
    s1 = state & ~(1 << q)                # remove particle at q
    if (s1 >> p) & 1:                     # c^dag_p needs p empty
        return None, 0.0
    sgn *= jw_sign(s1, p)
    s2 = s1 | (1 << p)                    # add particle at p
    return s2, sgn


def build_hop_matrix(sector_states, index_of, p, q):
    """Matrix of c^dag_p c_q on the sector (non-Hermitian on its own)."""
    dim = len(sector_states)
    M = np.zeros((dim, dim), dtype=complex)
    for i, s in enumerate(sector_states):
        s2, amp = apply_cdag_c(s, p, q)
        if s2 is not None and s2 in index_of:
            M[index_of[s2], i] += amp
    return M


def build_coherence_operators(sector_states, index_of, mK, mKp):
    """
    Intervalley coherence operators as sector matrices:
        T_x = c^dag_{mK} c_{mKp} + h.c.
        T_y = -i ( c^dag_{mK} c_{mKp} - h.c. )
    (the x, y components of the valley pseudospin; T_z is valley_pol_matrix).
    """
    A = build_hop_matrix(sector_states, index_of, mK, mKp)   # c^dag_{mK} c_{mKp}
    Adag = A.conj().T
    Tx = A + Adag
    Ty = -1j * (A - Adag)
    return Tx, Ty


#Ground-space projection and the valley-symmetric state
def ground_space_polarization(H_sec, Pmat, deg_tol=1e-6):
    """
    Diagonalize the valley polarization T_z within the (degenerate) ground
    space of H_sec. Returns (E0, degeneracy, polarization eigenvalues).
    Basis-independent: the eigenvalues do not depend on how eigh orients the
    degenerate ground vectors.
    """
    evals, evecs = np.linalg.eigh(H_sec)
    E0 = evals[0]
    deg = int(np.sum(np.abs(evals - E0) < deg_tol))
    G = evecs[:, :deg]
    P2 = G.conj().T @ Pmat @ G
    return E0, deg, np.linalg.eigvalsh(P2)


def project_into_ground_space(H_sec, Op, deg_tol=1e-6):
    """
    Project an arbitrary sector operator Op into the ground space of H_sec and
    return its eigenvalues. Used for T_x, T_y, T_z alike.
    """
    evals, evecs = np.linalg.eigh(H_sec)
    E0 = evals[0]
    deg = int(np.sum(np.abs(evals - E0) < deg_tol))
    G = evecs[:, :deg]
    Op2 = G.conj().T @ Op @ G
    return np.linalg.eigvalsh(Op2), np.linalg.norm(Op2)


def symmetric_state_energy(H_sec, sector_states, index_of, f1, f2, mK, mKp):
    """
    Energy of the valley-symmetric state: the two lower modes (f1, f2) filled,
    the extra particle in the equal superposition of the two manifold modes
    (mK, mKp). Its gap above the ground energy is the order parameter that
    lifts off zero as valley polarization sets in.
    """
    sK  = (1 << f1) | (1 << f2) | (1 << mK)
    sKp = (1 << f1) | (1 << f2) | (1 << mKp)
    v = np.zeros(len(sector_states), dtype=complex)
    v[index_of[sK]]  = 1 / np.sqrt(2)
    v[index_of[sKp]] = 1 / np.sqrt(2)
    return np.real(v.conj() @ H_sec @ v)


def commutator_norm(A, B):
    """Max-abs element of [A, B]; ~0 means the operators commute."""
    return np.max(np.abs(A @ B - B @ A))

