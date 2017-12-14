import numpy as np
import neighbourlist


def LJ_potential_ij(r, sigmaii, epsilonii, sigmajj, epsilonjj):
    
    sigmaij = 0.5 * (sigmai + sigmaj) # Lorentz-Berthelot combining rules: https://en.wikipedia.org/wiki/Combining_rules
    epsilonij = np.sqrt(epsiloni * epsilonj)

    ###############
    ##### use neighbourlist instead!!!!!! #####
    r_cut = 2.5*sigmaij
    if r >= r_cut:
        return 0.0
    ###!!!!!!!!!!!!!!!!!!!!!!####################
    ###############

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

def LJ_potential(box):
    '''Computes the total Lennard Jones potential of the system configuration of *box*.
    
    arguments:
    	box (Box): 	a box object, which has the (relevant) attributes *box.particles*, a numpy
    				array of Particle objects, each with their own position and sigma/epsilon, 
    				*box.dimension*, an int giving the dimension of the box (1d-3d), and *box.size*, 
    				a *dimension*-dimensional numpy array of float which gives the size of the box in each direction 
    '''

    ############
    # neighbourlist stuff
    # ???
    ############

    LJpot = 0.0 
    for i, particlei in enumerate(box.particles):
        LJpot_i = 0.0
        for particlej in box.particles[i+1:]:
            r = np.linalg.norm(enforce_pbc(particlei.position - particlej.position, box.size))
            LJpot_i += LJ_potential_ij(r, particlei.sigma, particlei.epsilon, 
            							  particlej.sigma, particlej.epsilon)
        LJpot += LJpot_i
    return LJpot