import numpy as np
import neighbourlist
import system

def LJ_potential_ij(r, sigmaii, epsilonii, sigmajj, epsilonjj, r_c, r_s):
    
    sigmaij = 0.5 * (sigmaii + sigmajj) # Lorentz-Berthelot: https://en.wikipedia.org/wiki/Combining_rules
    epsilonij = np.sqrt(epsilonii * epsilonjj)

    r_cut = r_c + r_s

    q = (sigmaij / r)**6
    q_c = (sigmaij/r_cut)**6
    return (4.0 * epsilonij * q * (q - 1.0)  # return lennard jones interaction energy between particle i and particle j
           - 4.0 * epsilonij * q_c *(q_c - 1.0)) # subtract value of potential at r_cut to avoid discontinuity # precompute!!!!!!


def LJ_potential(box, r_c, r_s):    
    '''Computes the total Lennard Jones potential of the system configuration of *box*.
    
    arguments:
        box (Box): object of the class Box. Includes the above 2 variables in the object... 
                    (will be used when running the actual simulation, the 2 others (possibly) when testing)
        r_c (float): cutoff radius for LJ potential
        r_s (float): size of skin region for LJ potential
    '''

    if box.particles[0].neighbourlist is None:
        box.compute_LJneighbourlist(r_c, r_s)

    LJpot = 0.0
    for particlei in box.particles:
        for particlej in particlei.neighbourlist:
            r = np.linalg.norm(particlei.position - particlej.position)
            LJpot += LJ_potential_ij(r, particlei.sigmaLJ, particlei.epsilonLJ, 
                                        particlej.sigmaLJ, particlej.epsilonLJ, r_c, r_s)
    return LJpot/2 # fix? find cleverer solution