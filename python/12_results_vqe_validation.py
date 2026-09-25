#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import EfficientSU2
from qiskit.primitives import StatevectorEstimator
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA


def number_fermionic(n_modes):
    """Total particle-number operator as a FermionicOp."""
    data = {f"+_{m} -_{m}": 1.0 for m in range(n_modes)}
    return FermionicOp(data, num_spin_orbitals=n_modes)


def build_penalized_operator(qop, n_qubits, N_target, lam=5.0):
    """
    H_vqe = H + lam (N - N_target)^2, to confine VQE to the N_target sector.
    (N - N_target)^2 = N^2 - 2 N_target N + N_target^2 I.
    """
    mapper = JordanWignerMapper()
    N_op = mapper.map(number_fermionic(n_qubits))
    I_op = SparsePauliOp.from_list([("I" * n_qubits, 1.0)])
    penalty = lam * (N_op @ N_op - 2 * N_target * N_op + (N_target ** 2) * I_op)
    return (qop + penalty).simplify()


def run_vqe(H_vqe, n_qubits, reps=3, maxiter=800, seed=0):
    """Run EfficientSU2 + COBYLA VQE; return the optimizer result."""
    ansatz    = EfficientSU2(n_qubits, reps=reps, entanglement="full")
    estimator = StatevectorEstimator()
    optimizer = COBYLA(maxiter=maxiter)
    rng = np.random.default_rng(seed)
    x0 = rng.uniform(-np.pi, np.pi, ansatz.num_parameters)
    vqe = VQE(estimator, ansatz, optimizer, initial_point=x0)
    return vqe.compute_minimum_eigenvalue(H_vqe)


if __name__ == "__main__":
    U = 0.87
    n_qubits = N_modes

    # exact reference in the N=3 sector
    qop = to_qubit_op(H0_band, n, params, U)
    H_sec = sector_hamiltonian(qop, sector_states)
    exact = np.linalg.eigvalsh(H_sec)
    E_exact = exact[0]
    deg = int(np.sum(np.abs(exact - E_exact) < 1e-8))
    print(f"exact N=3 ground energy: {E_exact:.8f}   (degeneracy {deg})")

    # penalized VQE (confined to N=3)
    H_vqe = build_penalized_operator(qop, n_qubits, N_target, lam=5.0)
    result = run_vqe(H_vqe, n_qubits)
    E_vqe = result.eigenvalue.real
    print(f"VQE ground energy:       {E_vqe:.8f}")
    print(f"error |VQE - exact|:     {abs(E_vqe - E_exact):.2e}")

