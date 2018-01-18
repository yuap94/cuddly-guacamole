import numpy as np
import numpy.testing as npt
import lennardjones
import system
import pbc


def test_LJ_potential_ij(dim):

    # kb = 1.38064852*10**(-13) # N*Å/K (Boltzmann constant)

    sigma_argon = 3.405 # Å
    epsilon_argon = 119.8 # actually epsilon/kb (K)

    r = np.linalg.norm(np.random.rand(dim))
    r_c = 2.5*sigma_argon
    n_skip = 10
    width = r_c / 10
    r_skin = 2*n_skip*width

    r_cut = r_c + r_skin

    LJpot_ArAr_ref = (4*epsilon_argon*((sigma_argon/r)**(12)-(sigma_argon/r)**(6)) 
                    - 4*epsilon_argon*((sigma_argon/r_cut)**(12)-(sigma_argon/r_cut)**(6)))

    LJpot_ArAr = lennardjones.LJ_potential_ij(r, sigma_argon, epsilon_argon, sigma_argon, epsilon_argon, r_c, r_skin)
    npt.assert_allclose(LJpot_ArAr, LJpot_ArAr_ref, rtol = 1e-12)


def test_LJ_potential(dim):

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
    
    ourbox = system.Box(dimension = dim, size = boxsize, particles = particles, temp = 120.0)

    LJ_pot_ref = 0
    
    for particlej in [argon_2, argon_3, xenon_1, xenon_2, xenon_3]:
        r = np.linalg.norm(argon_1.position - particlej.position)
        if r < r_c + r_s:
            LJ_pot_ref += lennardjones.LJ_potential_ij(r, sigma_argon, epsilon_argon, 
                            particlej.sigmaLJ, particlej.epsilonLJ, r_c, r_s)

    for particlej in [argon_3, xenon_1, xenon_2, xenon_3]:
        r = np.linalg.norm(argon_2.position - particlej.position)
        if r < r_c + r_s:
            LJ_pot_ref += lennardjones.LJ_potential_ij(r, sigma_argon, epsilon_argon, 
                            particlej.sigmaLJ, particlej.epsilonLJ, r_c, r_s)

    for particlej in [xenon_1, xenon_2, xenon_3]:
        r = np.linalg.norm(argon_3.position - particlej.position)
        if r < r_c + r_s:
            LJ_pot_ref += lennardjones.LJ_potential_ij(r, sigma_argon, epsilon_argon, 
                            particlej.sigmaLJ, particlej.epsilonLJ, r_c, r_s)

    for particlej in [xenon_2, xenon_3]:
        r = np.linalg.norm(xenon_1.position - particlej.position)
        if r < r_c + r_s:
            LJ_pot_ref += lennardjones.LJ_potential_ij(r, sigma_xenon, epsilon_xenon, 
                            particlej.sigmaLJ, particlej.epsilonLJ, r_c, r_s)

    for particlej in [xenon_3]:
        r = np.linalg.norm(xenon_2.position - particlej.position)
        if r < r_c + r_s:
            LJ_pot_ref += lennardjones.LJ_potential_ij(r, sigma_xenon, epsilon_xenon, 
                            particlej.sigmaLJ, particlej.epsilonLJ, r_c, r_s)

    ourbox.compute_LJneighbourlist(r_c, r_s)
    ourbox.compute_LJ_potential(r_c, r_s)
    # npt.assert_almost_equal(LJ_pot_ref, ourbox.LJpotential)
    npt.assert_allclose(LJ_pot_ref, ourbox.LJpotential, rtol=1e-12)




