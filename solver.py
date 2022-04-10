import numpy as np
import matplotlib.pyplot as plt

#Load TOV class
from tov import *
from constants import *

#Here we load EoS to be used in calculations
eos = np.genfromtxt("eos_unified_IOPB.dat")

n_arr = eos[:,0]
e_arr = eos[:,1]
p_arr = eos[:,2]

#initialization of TOVsolver
tov_s = TOV(e_arr, p_arr)

m_arr = []
R_arr = []
C_arr = []
dens_arr = []

for dens_c in np.logspace(-0.0,3.7,300):
     R, M, prof = tov_s.solve(dens_c)
     m_arr.append(M)
     R_arr.append(R)
     C_arr.append(M*ModotTokm/R)
     dens_arr.append(dens_c)
     print(M, R, dens_c)

plt.plot(R_arr, m_arr)

plt.ylabel(r'${\rm M/M_\odot}$')
plt.xlabel(r'${\rm R~(km)}$')

plt.savefig('mr.png')
plt.show()
print(max(m_arr))
    
OUT = np.c_[R_arr, m_arr, dens_arr, C_arr];
np.savetxt('star_IOPB.dat', OUT, fmt='%1.6e')
