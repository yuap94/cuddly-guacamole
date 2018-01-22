#Mathematical derivation. with reference to  http://micro.stanford.edu/mediawiki/images/4/46/Ewald_notes.pdf 

import numpy as np
import neighbourlist
import system
from numpy.linalg import norm
from scipy.special import erfc
from scipy import exp, pi





def energy(position,  q, cell, .......):
    """
    Arguments:
    position : potential location
    q : list of charges........

    """
    Energy_short  = short_energy_sum(i, r, q, cell, alpha, cutoff_rspace)
    Energy_long   = long_energy_sum(i, r, q, invcell, alpha, cutoff_kspace, area)
    Energy_self   = self_energy_sum(i, r, q, cell, alpha)
    
    return Energy_short+Energy_long-Engery_self


#energy calculation formula. with reference to equation 39 in page 7 of the pdf(link found in first line of this file). 
def short_energy_sum (i, r, q, cell, alpha, cutoff_rspace):
    
    for j in range(0, len(q)):
        rij = r[i, :] - r[j, :]
   
        
        
    return Vr


def long_energy_sum(i, r, q, cell, alpha, cutoff_kspace, area):
    long_energy_sum = 0
    for j in range(0, len(q)):
        rij = r[i, :] - r[j, :]
       
   
            
    long_pre = 1/(2*V*epsilon_0)     #prefactor of the long-ranged term         
    k                                #???reciprocal vector need to be defined 
    k2 = k**2 
    sigma2 = sigma** 2 
    s_k = sum( q*np.exp(k*r)         #structure factor 
    s_k2= s_k**2                     #sqaure of structure factor    

    midpart = np.sum(exp(-sigma2*k2/2)/k2) 
     
        long_energy_sum = long_pre*midpart*s_k2
    return long_energy_sum


def self_energy_sum(i, r, q, cell, alpha):
   
    for j in range(0, len(q)):
        if i == j:
            
                 
    self_pre = 1/(4*pi*eplison_0*sigma*np.sqrt(2*np.pi))     #prefactor of the self-term
    q2 = q** 2                                                  #???particle charge to be defined 
              
        self_energy_sum = self_pre*np.sum(q2) 
         
    return self_energy_sum 
