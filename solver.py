import numpy as np
import matplotlib.pyplot as plt

#Load TOV class
from tov import *
from constants import *

#Here we load EoS to be used in calculations
eos = np.genfromtxt("eos_unified_TM1.dat")

n_arr = eos[:,0]
e_arr = eos[:,1]
p_arr = eos[:,2]

#initialization of TOVsolver
tov_s = TOV(e_arr, p_arr)

#m_arr = []
#R_arr = []
#P_arr = []

#for dens_c in np.logspace(-0.0,3.7,300):
#     R, M, prof = tov_s.solve(dens_c)
#     m_arr.append(M)
#     R_arr.append(R)
#     print(M, R, dens_c)

#plt.plot(R_arr, m_arr)

#plt.ylabel(r'${\rm M/M_\odot}$')
#plt.xlabel(r'${\rm R~(km)}$')

#plt.savefig('mr.png')
#plt.show()
#print(max(m_arr))


# For canonical star
#R, M, prof = tov_s.solve(365.5)
#plt.plot(prof[0], prof[1])
#plt.show()

#Here we calculate tidal propertis of the NS family

m_arr = []
R_arr = []
den_arr = []
C_arr = []
k2_arr = []
L_arr = []

for dens_c in np.logspace(-0.0,3.7,150):
    NS_prop = tov_s.solve_tidal(dens_c, rmax=25)
    R, M, C, k2, L = NS_prop[0], NS_prop[1], NS_prop[2], NS_prop[3], NS_prop[4]
    print(R, M, C, k2, L)

    m_arr.append(M)
    R_arr.append(R)
    den_arr.append(dens_c)
    C_arr.append(C)
    L_arr.append(L)
    k2_arr.append(k2)
    
OUT = np.c_[R_arr, m_arr, den_arr, C_arr, k2_arr, L_arr];
np.savetxt('star_TM1.dat', OUT, fmt='%1.6e')


#fig = plt.figure(figsize=(10,8))

#ax1 = fig.add_subplot(221)

#ax1.plot(R_arr, m_arr)
#ax1.set_xlabel(r'${\rm R~(km)}$')
#ax1.set_ylabel(r'${\rm M/M_\odot}$')

#ax2 = fig.add_subplot(222)

#ax2.plot(den_arr, m_arr)
#ax2.set_ylabel(r'${\rm M/M_\odot}$')
#ax2.set_xlabel(r'${\rm den~(MeV/fm^{-3})}$')

#ax3 = fig.add_subplot(223)

#ax3.plot(C_arr, k2_arr)
#ax3.set_xlabel(r'${\rm R~(km)}$')
#ax3.set_ylabel(r'${\rm k_2}$')

#ax4 = fig.add_subplot(224)

#ax4.plot(m_arr, L_arr)
#ax4.set_xlabel(r'${\rm M/M_\odot}$')
#ax4.set_ylabel(r'${\rm \Lambda}$')
#ax4.set_xlim(0,2.2)
#ax4.set_ylim(0,1000)

#plt.show()
