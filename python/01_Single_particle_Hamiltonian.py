#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np


r"""
For Bi2Te3, $\hbar v_F = 0.3$ eV nm and $a \approx 13$ nm, which gives $E_0 = \hbar v_F/a = 23.08$ meV. 
Here, all energies are scaled in units of $E_0$ and all lengths are scaled in units of a.
"""

def build_hamiltonian(n,k,params):

    r"""
    Construct the Hamiltonian for a single Dirac fermion
    in a C6-symmetric moiré lattice potential V(r) under a constant magnetization hz.

    Model:
        $H = (k × σ)_z + V(r) σ_0 + hz σ_3$

    where
        $V(r) = 2 V0 ∑_{j=1}^3 cos(G_j · r)$

    Parameters
    ----------
    n : int
        Plane-wave cutoff along both the reciprocal lattice directions.
    k : list of float
        Crystal momentum.
    params : list of float
        [hz,V] whe hz is uniform magnetization strength and V is lattice potential amplitude.

    Returns
    -------
    H : ndarray (complex)
        Hamiltonian matrix of dimension 2 * n^2. 
    """

    #Crystal momenta
    [kx,ky] = k;

    #Parameters
    [hz,V] = params

    #Reciprocal lattice vectors
    N = int((n+1)**2);

    R = []

    dR = [-2*np.pi,-2*np.pi/np.sqrt(3)]
    xylist = np.linspace(-n/2,n/2,n+1)
    for i in xylist:
        for j in xylist :
            value = [(i-j)*dR[0],(i+j)*dR[1]]
            R.append(value)

    #Constructing Hamiltonian
    H = np.zeros((2*N,2*N),complex)

    #Diagonal terms
    for l in range(n+1):
        for m in range(n+1):
            H[l*(n+1) + m, l*(n+1) + m] = hz
            H[l*(n+1)+m, l*(n+1)+m + N] = -(1.j*(kx-R[l*(n+1)+m][0])+ (ky - R[l*(n+1)+m][1]))

            H[l*(n+1) + m + N, l*(n+1) + m + N] = -hz
            H[l*(n+1)+ m + N, l*(n+1) + m ] = (1.j*(kx-R[l*(n+1)+m][0])- (ky - R[l*(n+1)+m][1]))

    # Hopping along first direction (x-direction)        
    for l in range(n):
        for m in range(n+1):
            H[l*(n+1) + m, (l+1)*(n+1) + m] = V
            H[(l+1)*(n+1) + m , l*(n+1) + m] = V 

            H[l*(n+1) + m + N, (l+1)*(n+1) + m + N] = V  
            H[(l+1)*(n+1) + m + N, l*(n+1) + m + N] = V 

    # Hopping along second direction (y-direction)
    for l in range(n+1):
        for m in range(n):
            H[l*(n+1) + m , l*(n+1) + m+1] = V 
            H[l*(n+1) + m+1, l*(n+1) + m] = V 
            H[l*(n+1) + m + N, l*(n+1) + m+1 + N] = V 
            H[l*(n+1) + m+1 + N, l*(n+1) + m + N] = V  

    # Hopping along third direction (diagonal)
    for l in range(n):
        for m in range(n):
            H[l*(n+1) + m, (l+1)*(n+1) + m+1] = V
            H[(l+1)*(n+1) + m+1, l*(n+1) + m] =  V

            H[l*(n+1) + m + N, (l+1)*(n+1) + m+1 + N] =  V
            H[(l+1)*(n+1) + m+1 + N, l*(n+1) + m + N] =  V 

    return H

#Find energy given a k value and a band index r
def get_energy(n,k,params,r):
    """
    Parameters
    ----------
    n : int
        Plane-wave cutoff along both the reciprocal lattice directions.
    k : list of float
        Crystal momentum.
    params : list of float
        [hz,V] whe hz is uniform magnetization strength and V is lattice potential amplitude.
    r: int
        band index
    """
    H = build_hamiltonian(n,k,params);
    return np.real(np.sort(np.linalg.eigvals(H))[r])

