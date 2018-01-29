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



def self_potential ##


#energy calculation formula. with reference to equation 39 in page 7 of the pdf(link found in first line of this file). 
def short_energy_sum (i, r, q, cell, alpha, cutoff_rspace):
    
    for j in range(0, len(q)):
        rij = r[i, :] - r[j, :]
   
        
        
    return Vr


def long_e (i, r, q, cell, alpha, cutoff_kspace, area):
    long_e = 0
    
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


def total_self_energy(i, r, q, cell, alpha):                 #general algo. done 
    self_energy = 0 
    self_potential_sum = 0 
    self_pre = 1/(4*pi*eplison_0*sigma*np.sqrt(2*np.pi))     #prefactor of the self-term
              
    for i in range(0, len(q)):                               # !! particle code to be corrected to system.  
        self_potential_sum += self_potential(i)              # !! def self_potential 
      
    total_self_energy = self_pre*self_potential_sum 
         
    return total_self_energy
