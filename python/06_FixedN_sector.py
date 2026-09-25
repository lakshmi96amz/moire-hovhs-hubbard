#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from itertools import combinations

#Non-interacting reference in the 3-particle sector
N_target = 3          # 3/8 filling : half filling case at E = EK
n_qubits = N_modes    # 8

def fixed_N_states(n_qubits, N_target):
    """Return sorted list of integer basis states with exactly N_target bits set."""
    states = []
    for occ in combinations(range(n_qubits), N_target):
        s = 0
        for m in occ:
            s |= (1 << m)
        states.append(s)
    return sorted(states)

sector_states = fixed_N_states(n_qubits, N_target)
dim_sector = len(sector_states)
index_of = {s: i for i, s in enumerate(sector_states)}
print(f"sector dimension C({n_qubits},{N_target}) = {dim_sector}")

def sector_hamiltonian(qubit_op, sector_states):
    """Dense Hamiltonian restricted to the fixed-N sector."""
    H_full = qubit_op.to_matrix()          # 256 x 256
    idx = np.array(sector_states, dtype=int)
    H_sec = H_full[np.ix_(idx, idx)]
    return H_sec

#Diagonalizing at U=0
qop0 = to_qubit_op(H0_band, n, params, U=0.0)
H_sec_U0 = sector_hamiltonian(qop0, sector_states)

# Hermitian by construction; use eigh for sorted real eigenvalues.
evals_U0,_  = np.linalg.eigh(H_sec_U0)

print("lowest 8 sector eigenvalues at U=0:")
print(np.round(evals_U0[:8], 8))

#Turning on interaction to U = 0.87
H_sec_U = sector_hamiltonian(qop, sector_states)
evals_U,_ = np.linalg.eigh(H_sec_U)

print("lowest 8 sector eigenvalues at U=0.87:")
print(np.round(evals_U[:8], 8))

