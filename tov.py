from scipy.interpolate import interp1d
from scipy.integrate import odeint
import numpy as np
from constants import *

class TOV:
  
  def __init__(self, en_arr, p_arr):

    en_arr *= MeVfm_3Tokm_2 # km^-2
    p_arr  *= MeVfm_3Tokm_2 # km^-2
    #cs_arr = np.gradient(p_arr, en_arr)

    sort_ind = np.argsort(p_arr)
    self.en_dens = interp1d(p_arr[sort_ind], en_arr[sort_ind], kind='linear')

    sort_ind = np.argsort(en_arr)
    self.press = interp1d(en_arr[sort_ind], p_arr[sort_ind], kind='linear')

    # sort_ind = np.argsort(p_arr)
    # self.cs2 = interp1d(p_arr[sort_ind], cs_arr[sort_ind], kind='linear')

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

    #dPhidr = 2 * (m + 4.0 * pi * r ** 3 * P)/ (r**2 - 2.0 * m * r)  # km^-3

    return [dPdr, dmdr] #, dPhidr]

  def dedp(self, r, R_dep):
    e_R, p_R, m_R = R_dep

    p = p_R(r)
    dp = p * 0.005

    el_3 = self.en_dens(p - 3 * dp)
    el_2 = self.en_dens(p - 2 * dp)
    el_1 = self.en_dens(p - 1 * dp)
    er_3 = self.en_dens(p + 3 * dp)
    er_2 = self.en_dens(p + 2 * dp)
    er_1 = self.en_dens(p + 1 * dp)
    de_dp = (-1 / 60 * el_3 + 3 / 20 * el_2 - 3 / 4 * el_1 + 3 / 4 * er_1 - 3 / 20 * er_2 + 1 / 60 * er_3) / dp

    return de_dp

  def love_eq(self, param, r, R_dep):

    beta, H = param  # beta in km, H in km^2
    e_R, p_R, m_R = R_dep  # e_R, p_R in km^-2, m_R in km

    try:
      dummy = p_R(r)
    except ValueError:
      return [100000, 100000]

    de_dp = self.dedp(r, R_dep)

    #cs_2 = self.cs2

    ell = (1.0 - 2.0 * m_R(r)/r)**-1  # no unit

    r2 = r**2
    r3 = r**3

    nu2 = (2.0 * (m_R(r) + 4.0 * pi * r3 * p_R(r))/ (r2 - 2.0 * m_R(r) * r))**2

    dbetadr = - beta * ( (2.0 / r) + ell * ((2.0 * m_R(r) / r2) + (4.0 * pi * r *(e_R(r) - p_R(r))))) \
              - H * (ell * ((4.0 * pi *(e_R(r) + p_R(r)) * de_dp) + (4.0 * pi * (5.0 * e_R(r))) \
              + (4.0 * pi *(9.0 * p_R(r))) - 6.0 / r2) - nu2)
    
    dHdr = beta
    
    return [dbetadr, dHdr]

  def solve(self, c_dens, rmax=30, rtol=1.0e-6, dmrel=1.0e-6, dr=1.0e-2):
        
    c_dens = c_dens*MeVfm_3Tokm_2    

    self.check_density(c_dens)

    r = np.arange(dr, rmax, dr) # km

    P = self.press(c_dens)  # km^-2

    eden = self.en_dens(P)  # km^-2

    m = 4.0 * pi * r[0] ** 3 * eden   # km

    psol = odeint(self.tov_eq, [P, m], r)

    p_R, m_R = psol[:,0], psol[:,1]

    # phi_R = 2.0*(m + 4.0 * pi * r ** 3 * p_R)/ (r**2 - 2.0 * m_R * r) * r
    # const = (0.5*np.log(1.0 - 2.0*m_R/r) - phi_R)
    # phi_R = phi_R + const

  
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

    # phi_R = phi_R[:ind]
    # Phi_R = phi_R[ind-1]
    
    e_R = self.en_dens(p_R) # km^-2
    
    return R, M*kmToModot, (r, m_R, e_R, p_R) # R in km, M in M_sun, r in km, m_R in km, and e_R, p_R in km^-2 

  def solve_tidal(self, c_dens, rmax=25, rtol=1.0e-4, dmrel=1.0e-6, dr=1.0e-2):
    
    R, M, R_dep  = self.solve(c_dens, rmax=rmax, rtol=rtol, dmrel=dmrel, dr=dr)
    # R in km, M in M_sun
    r, m_R, e_R, p_R  = R_dep  # r, m_R in km, and e_R, p_R in km^-2 

    M *= ModotTokm

    e_R = interp1d(r, e_R, kind='linear')
    p_R = interp1d(r, p_R, kind='linear')
    m_R = interp1d(r, m_R, kind='linear')
  
    beta0 = 2.0 * r[0] # km
    H0 = r[0] ** 2.0   # km^2

    solution = odeint(self.love_eq, [beta0, H0], r, args=([e_R, p_R, m_R],), rtol=rtol)

    beta = solution[-1, 0]
    H = solution[-1, 1]

    y = R * beta / H

    C = M / R

    k2 = 8. / 5. * C ** 5 * (1 - 2 * C) ** 2 * (2 + 2 * C * (y - 1) - y) * (
         2 * C * (6 - 3 * y + 3 * C * (5 * y - 8)) + 4 * C ** 3 * (
         13 - 11 * y + C * (3 * y - 2) + 2 * C ** 2 * (1 + y)) + 3 * (1 - 2 * C) ** 2 * (2 - y + 2 * C * (y - 1)) * (
         np.log(1 - 2 * C))) ** (-1)

    Lamb = ((2./3.)*k2)/C**5     

    return np.array([R, M*kmToModot, C, k2, Lamb]) #, y, beta, H])
