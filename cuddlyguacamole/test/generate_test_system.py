# not currently in use...
import numpy as np
import numpy.testing as npt
import system
import pbc


def gen_test_box(dim):
    # kb = 1.38064852*10**(-13) # N*Å/K (Boltzmann constant)

    sigma_argon = 3.405 # Å
    epsilon_argon = 119.8 # actually epsilon/kb (K)

    sigma_xenon = 4.07 # Å
    epsilon_xenon = 225.3 # actually epsilon/kb (K)

    boxsize = np.ones(dim)*np.maximum(sigma_argon, sigma_xenon)*20

    r_c = 2.5*0.5*(sigma_argon+sigma_xenon)
    n_skip = 5
    width = r_c / (n_skip*10)
    r_s = 2*n_skip*width

    # pos1 = pbc.enforce_pbc(np.random.randn(dim)/4, boxsize)
    # pos2 = pbc.enforce_pbc(np.random.randn(dim)/4, boxsize)
    # pos3 = pbc.enforce_pbc(np.random.randn(dim)/4, boxsize)
    # pos4 = pbc.enforce_pbc(np.random.randn(dim)/4, boxsize)
    # pos5 = pbc.enforce_pbc(np.random.randn(dim)/4, boxsize)
    # pos6 = pbc.enforce_pbc(np.random.randn(dim)/4, boxsize)

    pos1 = pbc.enforce_pbc(np.random.randn(dim), boxsize)
    pos2 = pbc.enforce_pbc(pos1 + r_c * np.random.randn(dim), boxsize)
    pos3 = pbc.enforce_pbc(pos2 + r_c * np.random.randn(dim), boxsize)
    pos4 = pbc.enforce_pbc(pos3 + r_c * np.random.randn(dim), boxsize)
    pos5 = pbc.enforce_pbc(pos4 + r_c * np.random.randn(dim), boxsize)
    pos6 = pbc.enforce_pbc(pos5 + r_c * np.random.randn(dim), boxsize)
    
    argon_1 = system.Particle(position = pos1, charge = 0, sigmaLJ = sigma_argon, epsilonLJ = epsilon_argon)
    argon_2 = system.Particle(position = pos2, charge = 0, sigmaLJ = sigma_argon, epsilonLJ = epsilon_argon)
    argon_3 = system.Particle(position = pos3, charge = 0, sigmaLJ = sigma_argon, epsilonLJ = epsilon_argon)
    xenon_1 = system.Particle(position = pos4, charge = 0, sigmaLJ = sigma_xenon, epsilonLJ = epsilon_xenon)
    xenon_2 = system.Particle(position = pos5, charge = 0, sigmaLJ = sigma_xenon, epsilonLJ = epsilon_xenon)
    xenon_3 = system.Particle(position = pos6, charge = 0, sigmaLJ = sigma_xenon, epsilonLJ = epsilon_xenon)

    particles = [argon_1, argon_2, argon_3, xenon_1, xenon_2, xenon_3]
    
    return system.Box(dimension = dim, size = boxsize, particles = particles, temp = 120.0)