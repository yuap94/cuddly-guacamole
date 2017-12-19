import numpy as np
import lennardjones
import energy

def mcmc_step(positions_curr, box, potLJ_curr=None, width=0.2, temp = 273.15, update_nblist=True):

    kb = 1.38064852*10**(-23) # boltzmann constant

    if pot_LJ_curr is None:
        pot_LJ_curr = LJ_potential(box = box)

    if update_nblist:
        ####
        ###############
        ##########
        #TO DO!!!!!!!!!
        ##############
        ###############
        ##########

    positions_trial = (positions_curr + width * 
                      np.random.randn(*np.asarray(positions_curr).shape)/4) #randn -> std norm. dist, divide by 4 to keep results mostly within (-0.5,0.5)
    potLJ_trial = LJ_potential(positions = positions_trial, boxsize = box.size)

    if np.random.rand() < min(1,np.exp(-(potLJ_trial - potLJ_curr)/kb*box.temp)):
        return positions_trial, potLJ_trial
    return positions_curr, potLJ_curr

def mcmc(box, n_steps, width=0.2, n_skip=1, n_reuse_nblist = 1, 
            save_system_history = False, r_c = 2.5, r_s = 0.1):
    '''Metropolis MCMC simulation of the movement of the particles within *box*.
    NB: Currently only using LJ potential.
    
    arguments:
        box (Box): our box object 
        n_steps (int): how many mcmc steps in total to make
        width (float): max (??) size of each mcmc step
        n_skip (int): we only save the system state for every n_skip steps
        n_reuse_nblist (int): how many mcmc steps we do before we the update the 
                                verlet neighbourlist of the particles in box
        save_system_history (bool): save the history of system states during the mcmc simulation? 
    '''

    # Store initial position for each particle in list
    positions_initial = [] 
    for particle in box.particles:
        positions_initial.append(particle.position)
    # positions_initial = box.positions #<----  fix  

    # Store initial total LJ potential:
    potLJ_initial = LJ_potential(box, r_c, r_s) 

    # Save intial system configuration (positions & potential):
    if save_system_history:
        positions_history = [positions_initial] 
        potLJ_history = [potLJ_initial] 


    positions_curr, potLJ_curr = positions_initial, potLJ_initial 
    for i in range(int(n_steps/n_skip)):
        for j in range(n_skip):
            if mod(i+j,n_reuse_nblist) == 0:
                update_nblist = True
            positions_curr, potLJ_curr = mcmc_step(positions_curr, box, potLJ_curr, width, box.temp, update_nblist)
            update_nblist = False


        if save_system_history:
            positions_history.append([positions_curr])
            potLJ_history.append([potLJ_curr])

    for i, particlei in enumerate(box.particles): # update to final system config before returning box
        particlei.position = pos_curr[i] # or box.particles[i].position = pos_curr[i] ??? how do pointers work exactly in python?

    return box, positions_history, potLJ_history








