#Mathematical derivation. with reference to  http://micro.stanford.edu/mediawiki/images/4/46/Ewald_notes.pdf 

import numpy as np
import neighbourlist
import system
from numpy.linalg import norm
from scipy.special import erfc
from scipy import exp, pi


def energy(position,  q, cell, .......):

    Energy_short  = short_energy(i, r, q, cell, cutoff_rspace)m
    Energy_long   = long_energy(i, r, q, cutoff_kspace,)
    Energy_self   = self_energy(i, r, q, cell)
    
    return Energy_short+Energy_long-Engery_self




def self_potential ##

#energy calculation formula. with reference to equation 39 in page 7 of the pdf(link found in first line of this file). 
def short_energy_sum (i, r, q, cell, alpha, cutoff_rspace):
    
    for j in range(0, len(q)):
        
   
        
        
    return Vr


def k_energy (i, r, q, cutoff_kspace):
    #prefactor        
    pre_fac = 1/(2*V*epsilon_0)  
    
    result = np.zeros (n) 
    
    for i in range 

    
    
    
    #Reciprocal vector 
    k = 2.0*np.pi*np.array(, , )                               
    k2 = k**2 
    
    #Sigma = r_c/ sqrt(2*p)
    sigma2 = sigma** 2
    
    #Structure factor
    s_k = sum( q*np.exp(k*r))          
    s_k2= s_k**2              

    e_part = np.sum(exp(-sigma2*k2/2)/k2) 
     
        l_r_e = pre_fac*s_k2*e_part
        
    return l_r_e

def str_fac (self) :
    s_k=0 
    for i in range ():
        s_k += charge(i)*exp(k*position(i))
    return s_k 


















def total_s_energy(self):                                 #general algo. done 
    s_energy = 0 
    s_total_potential = 0 
    s_pre = 1/(4*pi*epsilon_0*sigma*np.sqrt(2*np.pi))     #prefactor of the self-term
              
    for i in range(0, len(q)):                            # ==1  particle code to be corrected to system.  
        total_s_potential += s_potential(i)               # !! def self_potential 
      
    total_s_energy = s_pre*s_total_potential
         
    return total_s_energy


