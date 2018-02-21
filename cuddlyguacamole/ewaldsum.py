#Mathematical derivation. with reference to  http://micro.stanford.edu/mediawiki/images/4/46/Ewald_notes.pdf 

import numpy as np
import neighbourlist
import system
from numpy.linalg import norm
from scipy.special import erfc
from scipy import exp, pi


def energy(position,  q, cell, .......):

    Energy_short  = short_energy(i, r, q, cell, cutoff_rspace)m
    Energy_long   = long_energy(i, r, q, cutoff_kspace)
    Energy_self   = self_energy(i, r, q, cell)
    
    return Energy_short+Energy_long-Engery_self




def self_potential ##

#energy calculation formula. with reference to equation 39 in page 7 of the pdf(link found in first line of this file). 
def short_energy_sum (i, r, q, cell, alpha, cutoff_rspace):
    
    for j in range(0, len(q)):
        
   
        
        
    return Vr






def k_energy (i, r, q, cutoff_kspace):
    
    
    
    for k_i in range (-k_c,k_c+1)
        for k_j in range (-k_c,k_c+1)
            for k_k in range (-k_c,k_c+1) 
                #Reciprocal vector 
                k = 2.0*np.pi*np.array([k_i / (box[0]), k_j / (box[1]), k_k / (box[2])] )                               
                k2 = k*k  
    
    #prefactor, V value ???
    pre_fac = 1/(2*V*epsilon_0)  

    #Structure factor
    
    s_k = sum( q*np.exp(k*r))          
    s_k2= s_k**2              

    exp_term = np.sum(exp(-sigma2*k2/2)/k2) 
     
         k_space_energy= pre_fac*s_k2*exp_term
        
    return k_space_energy


def k_cut_off (self) :
    #accuracy set to 1e-6 gives p 
    p = ????
    k_c = 2*p/r_cut 
    #Sigma
    sigma = r_c / np.sqrt (2*p) 
    sigma2 = sigma*sigma 
    

                      
 
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


