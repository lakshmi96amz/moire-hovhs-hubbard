#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#U=0 case
H0_band = build_H0_band(n, N_b, params);

#Validation

#Hermiticity of the assembled operator at U != 0
qop = to_qubit_op(H0_band, n, params, U=0.87)
diff = (qop - qop.adjoint()).simplify()
herm = np.max(np.abs(diff.coeffs)) if len(diff.coeffs) else 0.0
print(f"assembled non-Hermiticity: {herm:.2e}")  

#coefficients real
print(f"max imag coeff: {np.max(np.abs(qop.coeffs.imag)):.2e}") 

#U=0 recovers the non-interacting operator
qop0 = to_qubit_op(H0_band, n, params, U=0.0)
d0   = (qop0 - qop0.adjoint()).simplify()
print(f"U=0 non-Hermiticity: {np.max(np.abs(d0.coeffs)) if len(d0.coeffs) else 0.0:.2e}")

#positive-semidefinite interaction: <psi|H_int|psi> >= 0 for a repulsive Hubbard U
Hint_only = JordanWignerMapper().map(build_H_int(n, params, U=0.87))
Hmat = Hint_only.to_matrix()
rng  = np.random.default_rng(0)
vals = []
for _ in range(20):
    psi = rng.standard_normal(Hmat.shape[0]) + 1j*rng.standard_normal(Hmat.shape[0])
    psi /= np.linalg.norm(psi)
    vals.append((psi.conj() @ Hmat @ psi).real)
print(f"min <psi|H_int|psi>: {min(vals):.3e}")


