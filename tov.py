from scipy.interpolate import interp1d
from scipy.integrate import odeint
import numpy as np
from constants import *

class TOV:
  
  def __init__(self, en_arr, p_arr):

    en_arr *= MeVfm_3Tokm_2 # km^-2
    p_arr  *= MeVfm_3Tokm_2 # km^-2

    sort_ind = np.argsort(p_arr)
    self.en_dens = interp1d(p_arr[sort_ind], en_arr[sort_ind], kind='linear')

    sort_ind = np.argsort(en_arr)
    self.press = interp1d(en_arr[sort_ind], p_arr[sort_ind], kind='linear')

    self.min_dens = np.min(en_arr)
    self.max_dens = np.max(en_arr)

    self.min_p = np.min(p_arr)
    self.max_p = np.max(p_arr)
    
  def check_density(self, dens):
    if dens < self.min_dens or dens > self.max_dens:
      raise Exception('Central density: %8.4E is outside of the EoS range. \n' 
                        %(dens) +
                        'min density is: %8.4E, max density is:%8.4E'
                         %(self.min_dens, self.max_dens))  

  def tov_eq(self, y, r):
    
    P, m = y   # Phi in no unit, P in km^-2, m in km and r in km

    if P < self.min_p or P > self.max_p:
        
      return [0., 0.]

    eden = self.en_dens(P)  # km^-2

    dPdr = -(eden + P) * (m + 4.0 * pi * r ** 3 * P)/ (r**2 - 2.0 * m * r) # km^-3

    dmdr = 4.0 * pi * r ** 2 * eden  # no unit

    return [dPdr, dmdr]

  def solve(self, c_dens, rmax=30, rtol=1.0e-6, dmrel=1.0e-6, dr=1.0e-2):
        
    c_dens = c_dens*MeVfm_3Tokm_2    

    self.check_density(c_dens)

    r = np.arange(dr, rmax, dr) # km

    P = self.press(c_dens)  # km^-2

    eden = self.en_dens(P)  # km^-2

    m = 4.0 * pi * r[0] ** 3 * eden   # km

    psol = odeint(self.tov_eq, [P, m], r)

    p_R, m_R = psol[:,0], psol[:,1]

    diff = (m_R[1:] - m_R[:-1])/m_R[1:]
    ind = -1
    for i, dm in enumerate(diff):
      if dm < dmrel and m_R[i] != 0:
        ind = i
        break

    M = m_R[ind - 1]
    R = r[ind - 1]
    
    r   = r[:ind]   # km 
    p_R = p_R[:ind] # km^-2 
    m_R = m_R[:ind] # km
    
    e_R = self.en_dens(p_R) # km^-2
    
    return R, M*kmToModot, (r, m_R, e_R, p_R) # R in km, M in M_sun, r in km, m_R in km, and e_R, p_R in km^-2 
