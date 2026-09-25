#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np

#Non-interacting spectrum and manifold identification
#
# Inspect the single-particle spectrum to locate the HOVHS manifold and confirm
# the Fermi level sits at E_K. This defines the half-filling target (N=3) and
# the filled/manifold mode structure used by all downstream analysis.

H0_band = build_H0_band(n, N_b, params);
eps = np.sort(np.diag(H0_band).real)
EK  = get_energy(n, K, params, 121)      # HOVHS energy (manifold band at K)

print("single-particle energies (sorted):", np.round(eps, 6))
print("E_K =", round(EK, 6))

n_below = int(np.sum(eps < EK - 1e-6))
n_at    = int(np.sum(np.abs(eps - EK) < 1e-6))
print("modes below E_K:", n_below)
print("modes at E_K (manifold degeneracy):", n_at)

# manifold-to-nearest-band gap Delta (sets the clean-regime ceiling for U):
# distance from the manifold energy E_K to the nearest single-particle level
# outside the manifold, on whichever side is closer.
below = eps[eps < EK - 1e-6]
above = eps[eps > EK + 1e-6]
gap_below = EK - below.max() if len(below) else np.inf
gap_above = above.min() - EK if len(above) else np.inf
Delta = min(gap_below, gap_above)
print("manifold-to-nearest-band gap Delta =", round(Delta, 6))

