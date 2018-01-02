import numpy as np
import lennardjones
import neighbourlist
#import energy # for coloumb energy
import system
import numpy.testing as npt
import copy
import warnings


def mcmc_step(box, width, r_cut, r_skin, update_nblist):

    # kb = 1.38064852*10**(-13) # N*Å/K (Boltzmann constant)

    box_trial = copy.deepcopy(box)
    trial_step = width * np.random.randn(*np.asarray(box_trial.positions).shape)/4 #randn -> std norm. dist, divide by 4 to keep results mostly within (-0.5, 0.5)

    for i, particle in enumerate(box_trial.particles):
        box_trial.particles[i].position = system.enforce_pbc(box_trial.particles[i].position + trial_step[i], box_trial.size) 
        # print(box_trial.particles[0].position)

    box_trial.make_positions_list()
    box_trial.compute_LJ_potential(r_cut, r_skin)
    if update_nblist:
        box_trial.compute_LJneighbourlist(r_cut, r_skin) 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        acceptance_prob = min(1.0, np.exp(-(box_trial.LJpotential - box.LJpotential)/box.temp)) # note: potential is actually potential/kb

    # if (box_trial.LJpotential>box.LJpotential):
    #     print("increase = " + repr(box_trial.LJpotential-box.LJpotential))
    #     print("acceptance prob. = " + repr(acceptance_prob))

    # print(acceptance_prob)
    if np.random.rand() < acceptance_prob:
        return box_trial, trial_step, acceptance_prob # return trial_step and acceptance_prob to use in unit testing
    return box, trial_step, acceptance_prob


def mcmc(box, n_steps, width, n_skip, n_reuse_nblist, 
            save_system_history, r_cut_LJ, r_skin_LJ):
    '''Metropolis MCMC simulation of the movement of the particles within *box*.
    NB: Currently only using LJ potential. Returns the box object with new system configuration
    and a system configuration history.
    
    arguments:
        box (Box): our box object 
        n_steps (int): how many mcmc steps in total to make
        width (float): approximate "max" size of each mcmc step
        n_skip (int): we only save the system state for every n_skip steps
        n_reuse_nblist (int): how many mcmc steps to do before updating the LJ
                                verlet neighbourlist of the particles in *box*
        save_system_history (bool): save the history of system states during the mcmc simulation? 
        r_cut_LJ (float): cut-off radius for LJ neighbourlist computation
        r_skin_LJ (float): size of skin-reigon for LJ neighbourlist computation
    '''
    # Store initial position for each particle in list
    positions_history = [box.positions]             
    # (Compute then) Store initial total LJ potential:
    if box.LJpotential is None:
        box.compute_LJ_potential(r_cut_LJ, r_skin_LJ)
    potLJ_history = [box.LJpotential]  

    p_acc_vec = []
    for i in range(int(np.ceil(n_steps/n_skip))):
        for j in range(n_skip):
            # box_old = copy.deepcopy(box)
            if np.mod(i*n_skip+j, n_reuse_nblist+1) == 0:
                box, _, p_acc = mcmc_step(box, width, r_cut_LJ, r_skin_LJ, True) # mcmc acceptance prob p_acc used in testing
            else:
                box, _, p_acc = mcmc_step(box, width, r_cut_LJ, r_skin_LJ, False)
        # if (box.LJpotential>box_old.LJpotential):
        #     print("increase = " + repr(box.LJpotential-box_old.LJpotential))
        #     print("acceptance prob. = " + repr(p_acc))
        if save_system_history:
            positions_history.append(box.positions)
            potLJ_history.append(box.LJpotential)
            p_acc_vec.append(p_acc)

    return box, positions_history, potLJ_history, p_acc_vec # return p_acc_vec for use in testing



