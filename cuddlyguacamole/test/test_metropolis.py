import numpy as np
import numpy.testing as npt
import metropolis
import system
import test.generate_test_system
import lennardjones
import copy
import warnings
import pbc

def test_mcmc_step(dim, update_nblist):

    # kb = 1.38064852*10**(-13) # N*Å/K (Boltzmann constant)

    sigma_argon = 3.405 # Å
    epsilon_argon = 119.8 # # actually epsilon/kb (K)

    sigma_xenon = 4.07 # Å
    epsilon_xenon = 225.3 # actually epsilon/kb (K)

    boxsize = np.ones(dim)*np.maximum(sigma_argon, sigma_xenon)*20 # Å

    r_c = 2.5*0.5*(sigma_argon+sigma_xenon)
    n_skip = 5
    width = r_c / (n_skip*10)
    r_s = 2*n_skip*width

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

    temp = 120.0
    ourbox = system.Box(dimension = dim, size = boxsize, particles = particles, temp = temp)
    ourbox.compute_LJneighbourlist(r_c, r_s)
    ourbox.compute_LJ_potential(r_c, r_s)
    ourbox.make_positions_list()

    ourbox_trial = copy.deepcopy(ourbox)

    ourbox_post_mcmc_step, trial_step, _ = metropolis.mcmc_step(ourbox, width, r_c, r_s, update_nblist)

    for i, particle in enumerate(ourbox_trial.particles):
        ourbox_trial.particles[i].position = pbc.enforce_pbc(ourbox_trial.particles[i].position + trial_step[i], ourbox_trial.size)

    ourbox_trial.make_positions_list()
    ourbox_trial.compute_LJ_potential(r_c, r_s)
    if update_nblist:
        ourbox_trial.compute_LJneighbourlist(r_c, r_s)
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        prob_accept = np.minimum(1, np.exp(-(ourbox_trial.LJpotential - ourbox.LJpotential)/temp))
    
    if np.random.rand() < prob_accept:
        return ourbox_trial, ourbox_post_mcmc_step
    else:
        return ourbox, ourbox_post_mcmc_step





def test_mcmc(dim, n_steps, n_skip, n_reuse_nblist):

    # kb = 1.38064852*10**(-13) # N*Å/K (Boltzmann constant)
    sigma_argon = 3.405 # Å
    epsilon_argon = 119.8 # # actually epsilon/kb (K)
    sigma_xenon = 4.07 # Å
    epsilon_xenon = 225.3 # actually epsilon/kb (K)

    boxsize = np.ones(dim)*np.maximum(sigma_argon, sigma_xenon)*20 # Å

    r_c = 2.5*0.5*(sigma_argon+sigma_xenon)
    # n_skip = int(n_steps/10)
    width = r_c / (n_skip*10)
    r_s = 2*n_skip*width

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

    temp = 120.0
    ourbox = system.Box(dimension = dim, size = boxsize, particles = particles, temp = temp)
    ourbox.compute_LJneighbourlist(r_c, r_s)
    ourbox.compute_LJ_potential(r_c, r_s)
    ourbox.make_positions_list()

    ourbox_post_mcmc, positions_history, potLJ_history, p_acc_vec = metropolis.mcmc(ourbox, n_steps, width, n_skip, n_reuse_nblist, True, r_c, r_s)
   
    for positions in positions_history:
        temp_box = system.Box(dimension = dim, size = boxsize, particles = particles, temp = temp)
        for i in range(len(positions_history)):
            for j in range(len(temp_box.particles)):
                temp_box.particles[j].position = positions_history[i][j]
            temp_box.compute_LJ_potential(r_c, r_s)
            npt.assert_almost_equal(temp_box.LJpotential, potLJ_history[i], decimal = 3)
    # print(potLJ_history)









