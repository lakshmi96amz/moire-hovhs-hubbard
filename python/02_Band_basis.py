#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from itertools import product


#Gets U vectors for a given k value, params and sigma
def get_U(n, params, kvec, sigma):
    r"""
    n       : int -- Plane-wave cutoff along both the reciprocal lattice directions
    params  : list of float -- [hz,V] whe hz is uniform magnetization strength and V is lattice potential amplitude.
    kvec    : list of float -- crystal momentum
    sigma   : int -- 0/1 indicates spin up/down

    Returns : U^{\sigma}(kvec)
    """
    H = build_hamiltonian(n,kvec,params);
    val,vec = np.linalg.eigh(H);
    idx = np.argsort(val);
    sorted_vec = vec[:,idx]

    if sigma == 0:
        return sorted_vec[0:N_vec,bands]
    else:
        return sorted_vec[N_vec : 2*N_vec,bands]

#Generating the overlap value M(k1,k2,r1,r2,sigma) = \sum_i (U^{\sigma}_{i,r1}(k1))^* U^{\sigma}_{i,r2}(k2)
def get_overlap(n, params, r1,r2,kvec1, kvec2,sigma):
    r"""
    n            : int -- Plane-wave cutoff along both the reciprocal lattice directions
    params       : list of float -- [hz,V] whe hz is uniform magnetization strength and V is lattice potential amplitude.
    r1,r2        : int -- band indices ranging from 0-3
    kvec1, kvec2 : list of floats -- crystal momenta
    sigma        : int -- 0/1 indicates spin up/down

    Returns      : \sum_i (U^{\sigma}_{i,r1}(k1))^* U^{\sigma}_{i,r2}(k2)
    """
    #Computng the corresponding eigenvectors for the overlap
    vec1 = get_U(n, params, kvec1, sigma)[:,r1];
    vec2 = get_U(n, params, kvec2, sigma)[:,r2];
    return np.dot(vec1.conj(), vec2)

# Primitive moiré reciprocal lattice vectors, read off from build_hamiltonian's R construction.
# R value = [(i-j)*dR0, (i+j)*dR1], so stepping i by 1 gives b1, stepping j by 1 gives b2.
dR0 = -2*np.pi
dR1 = -2*np.pi/np.sqrt(3)
b1 = np.array([ dR0,  dR1])   # (i -> i+1): (+1)*dR0 in x, (+1)*dR1 in y
b2 = np.array([-dR0,  dR1])   # (j -> j+1): (-1)*dR0 in x, (+1)*dR1 in y

def find_patch_index_folded(kvec4, Kpatch, b1, b2, nmax=3, atol=1e-6):
    """
    Return the patch index p such that kvec4 == Kpatch[p] + G
    for some moiré reciprocal lattice vector G = m1*b1 + m2*b2,
    with m1, m2 in [-nmax, nmax]. Returns None if no match.
    """
    kv = np.asarray(kvec4, dtype=float)
    for p, kp in enumerate(Kpatch):
        kp = np.asarray(kp, dtype=float)
        diff = kv - kp
        # Solve diff = m1*b1 + m2*b2 for integer m1, m2.
        # Build the 2x2 matrix [b1 b2] and invert.
        B = np.column_stack([b1, b2])
        m = np.linalg.solve(B, diff)
        # Check m is (near) integer within the search range.
        m_round = np.round(m)
        if (np.all(np.abs(m - m_round) < atol) and
                np.all(np.abs(m_round) <= nmax)):
            return p
    return None



