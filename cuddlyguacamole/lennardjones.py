import numpy as np
#import system?


def LJ_potential_ij(r, sigmaii, epsilonii, sigmajj, epsilonjj, r_c, r_s):
    
    sigmaij = 0.5 * (sigmai + sigmaj) # Lorentz-Berthelot: https://en.wikipedia.org/wiki/Combining_rules
    epsilonij = np.sqrt(epsiloni * epsilonj)

    r_cut = r_c + r_s

    q = (sigmaij / r)**6
    q_c = (sigmaij/r_cut)**6
    return (4.0 * epsilonij * q * (q - 1.0) 
    	   - 4.0 * epsilonij * q_c *(q_c - 1.0)) # subtract value of potential at r_cut to avoid discontinuity


def enforce_pbc(r_vec, boxsize):
    for i, length in enumerate(boxsize):
        while r_vec[i] >= 0.5 * length:
            r_vec[i] -= length
        while r_vec[i] < -0.5 * length:
            r_vec[i] += length
    return r_vec


def LJ_potential(box, r_c, r_s):	
    '''Computes the total Lennard Jones potential of the system configuration of *box*.
    
    arguments:
		#positions (python list of d-dimensional numpy arrays (d=1-3)): the i-th array in the
					list gives the position of the i-th particle in the system
		#boxsize (d-dimensional numpy array of float): gives the size of the box in each direction 
		box (Box): object of the class Box. Includes the above 2 variables in the object... 
					(will be used when running the actual simulation, the 2 others (possibly) when testing)
		r_c (float): cutoff radius for LJ potential
		r_s (float): size of skin region for LJ potential
    '''
    LJpot = 0.0 
    for particlei in box.particles:
        LJpot_i = 0.0
        for particlej in particlei.neighbourlist:
            r = np.linalg.norm(enforce_pbc(particlei.position - particlej.position, box.size))
            LJpot_i += LJ_potential_ij(r, particlei.sigma, particlei.epsilon, 
        								particlej.sigma, particlej.epsilon, r_c, r_s)
        LJpot += LJpot_i

    return LJpot/2