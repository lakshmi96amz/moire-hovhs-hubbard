#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Momentum cutoff
n = 10;

#Moir\'e potential and Zeeman field strengths
params = [0.5,1.38]; #[V,hz]

#Working with just 2 patch points: K and K' for the inter-valley momentum scattering.
K = [2*np.pi/3 , 2*np.pi/np.sqrt(3)];
Kp = [-2*np.pi/3 , -2*np.pi/np.sqrt(3)];
Kpatch = [K, Kp];
N_k = len(Kpatch);

#Chosen band indices for n = 10
bands = [120,121,122,123];
N_b = len(bands);

N_vec = (n+1)**2 ; #length of each eigenvector spin-wise

N_modes = N_k * N_b #No. of modes

