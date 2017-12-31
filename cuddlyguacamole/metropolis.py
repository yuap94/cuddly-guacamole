import numpy as np
import lennardjones
import neighbourlist
#import energy # for coloumb energy
import system

def mcmc_step(box, width, r_cut, r_skin, update_nblist):

    kb = 1.38064852*10**(-23) # boltzmann constant

    box_trial = box
    box_trial.positions = box_trial.positions + width * np.random.randn(*np.asarray(box_trial.positions).shape)/4
    for i in range(len(box_trial.positions)):
        box_trial.positions[i] = system.enforce_pbc(box_trial.positions[i], box_trial.size)
                        #randn -> std norm. dist, divide by 4 to keep results mostly within (-0.5, 0.5)
    if update_nblist:
        box.neighourlist = neighbourlist.verlet_neighbourlist(box, r_cut, r_skin)
        box_trial.neighourlist = neighbourlist.verlet_neighbourlist(box_trial, r_cut, r_skin) # fix? indeed necessary to also compute nblist for trial box?
    
    box_trial.LJpotential = lennardjones.LJ_potential(box_trial, r_cut, r_skin) 

    if np.random.rand() < min(1, np.exp(-(box_trial.LJpotential - box.LJpotential)/(kb*box.temp))):
        # print("hi")
        return box_trial
    return box


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
        r_sking_LJ (float): size of skin-reigon for LJ neighbourlist computation
    '''

    # Store initial position for each particle in list
    if box.positions is None:
        box.make_positions_list()

    # Store initial total LJ potential:
    if box.LJpotential is None:
        box.compute_LJ_potential(r_cut_LJ, r_skin_LJ)

    # Save intial system configuration (positions & potential):
    positions_history = [box.positions] 
    potLJ_history = [box.LJpotential] 
 
    for i in range(int(n_steps/n_skip)): # fix?!!!
        for j in range(n_skip):
            if np.mod(i+j, n_reuse_nblist) == 0:
                box = mcmc_step(box, width, r_cut_LJ, r_skin_LJ, True)
                for k in range(len(box.particles)): #update position for each particle object in box
                    box.particles[k].position = box.positions[k]
            else:
                box = mcmc_step(box, width, r_cut_LJ, r_skin_LJ, False)
                for k in range(len(box.particles)): #update position for each particle object in box
                    box.particles[k].position = box.positions[k]

        if save_system_history:
            positions_history.append(box.positions)
            potLJ_history.append(box.LJpotential)

    for i in range(len(box.particles)): #update position for each particle object in box before returning box
        box.particles[i].position = box.positions[i]

    return box, positions_history, potLJ_history



