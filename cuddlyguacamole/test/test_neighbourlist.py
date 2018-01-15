import numpy as np
import numpy.testing as npt
import neighbourlist
import system    
import matplotlib.pyplot as plt
import pbc


def test_verlet_neighbourlist(dim):

    # kb = 1.38064852*10**(-13) # N*Å/K (Boltzmann constant)

    sigma_argon = 3.405 # Å
    epsilon_argon = 119.8 # actually epsilon/kb (K)

    sigma_xenon = 4.07 # Å
    epsilon_xenon = 225.3 # actually epsilon/kb (K)

    r_c = 2.5*0.5*(sigma_argon+sigma_xenon)
    width = r_c / 10
    n_skip = 10
    r_s = 2*n_skip*width

    boxsize = np.ones(dim)*np.maximum(sigma_argon, sigma_xenon)*20 # Å

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


    argon_1.neighbourlist = []
    for particlej in [argon_2, argon_3, xenon_1, xenon_2, xenon_3]:
        if np.linalg.norm(argon_1.position - particlej.position) < r_c + r_s:
            argon_1.neighbourlist.append(particlej)

    argon_2.neighbourlist = []
    for particlej in [argon_1, argon_3, xenon_1, xenon_2, xenon_3]:
        if np.linalg.norm(argon_2.position - particlej.position) < r_c + r_s:
            argon_2.neighbourlist.append(particlej)

    argon_3.neighbourlist = []
    for particlej in [argon_1, argon_2, xenon_1, xenon_2, xenon_3]:
        if np.linalg.norm(argon_3.position - particlej.position) < r_c + r_s:
            argon_3.neighbourlist.append(particlej)

    xenon_1.neighbourlist = []
    for particlej in [argon_1, argon_2, argon_3, xenon_2, xenon_3]:
        if np.linalg.norm(xenon_1.position - particlej.position) < r_c + r_s:
            xenon_1.neighbourlist.append(particlej)

    xenon_2.neighbourlist = []
    for particlej in [argon_1, argon_2, argon_3, xenon_1, xenon_3]:
        if np.linalg.norm(xenon_2.position - particlej.position) < r_c + r_s:
            xenon_2.neighbourlist.append(particlej)

    xenon_3.neighbourlist = []
    for particlej in [argon_1, argon_2, argon_3, xenon_1, xenon_2]:
        if np.linalg.norm(xenon_3.position - particlej.position) < r_c + r_s:
            xenon_3.neighbourlist.append(particlej)

    ourbox.compute_LJneighbourlist(r_c, r_s)

    for i in range(len(particles)):
        if particles[i].neighbourlist != ourbox.particles[i].neighbourlist:
            print("neighbourlists different")
        # print("particles[" + repr(i) + "] = \n" + repr(particles[i].neighbourlist))
        # print("ourbox.particles[" + repr(i) + "] = \n" + repr(ourbox.particles[i].neighbourlist))



# Neighbourlist/cut-off radius illustration:

def plot_neighbourlist():
    sigma_argon = 3.405 # Å
    epsilon_argon = 119.8 # actually epsilon/kb (K)

    sigma_xenon = 4.07 # Å
    epsilon_xenon = 225.3 # actually epsilon/kb (K)

    r_c = 2.5*0.5*(sigma_argon+sigma_xenon)
    width = r_c / 10
    n_skip = 10
    r_s = 2*n_skip*width

    boxsize = np.ones(2)*np.maximum(sigma_argon, sigma_xenon)*20 # Å

    pos1 = pbc.enforce_pbc(np.random.randn(2), boxsize)
    pos2 = pbc.enforce_pbc(pos1 + r_c * np.random.randn(2), boxsize)
    pos3 = pbc.enforce_pbc(pos2 + r_c * np.random.randn(2), boxsize)
    pos4 = pbc.enforce_pbc(pos3 + r_c * np.random.randn(2), boxsize)
    pos5 = pbc.enforce_pbc(pos4 + r_c * np.random.randn(2), boxsize)
    pos6 = pbc.enforce_pbc(pos5 + r_c * np.random.randn(2), boxsize)

    argon_1 = system.Particle(position = pos1, charge = 0, sigmaLJ = sigma_argon, epsilonLJ = epsilon_argon)
    argon_2 = system.Particle(position = pos2, charge = 0, sigmaLJ = sigma_argon, epsilonLJ = epsilon_argon)
    argon_3 = system.Particle(position = pos3, charge = 0, sigmaLJ = sigma_argon, epsilonLJ = epsilon_argon)
    xenon_1 = system.Particle(position = pos4, charge = 0, sigmaLJ = sigma_xenon, epsilonLJ = epsilon_xenon)
    xenon_2 = system.Particle(position = pos5, charge = 0, sigmaLJ = sigma_xenon, epsilonLJ = epsilon_xenon)
    xenon_3 = system.Particle(position = pos6, charge = 0, sigmaLJ = sigma_xenon, epsilonLJ = epsilon_xenon)

    particles = [argon_1, argon_2, argon_3, xenon_1, xenon_2, xenon_3]
    
    ourbox = system.Box(dimension = 2, size = boxsize, particles = particles, temp = 120.0)
    ourbox.compute_LJneighbourlist(r_c, r_s)


    j=1
    nbhood_circle = plt.Circle(ourbox.particles[j].position, r_c + r_s, color='b', fill=False)
    fig = plt.figure()

    ax = fig.add_subplot(111)
    ax.set_xlim((ourbox.center[0]-ourbox.size[0], ourbox.center[0]+ourbox.size[0]))
    ax.set_ylim((ourbox.center[1]-ourbox.size[1], ourbox.center[1]+ourbox.size[1]))

    ax.add_artist(nbhood_circle)

    nbpositionsj = [ourbox.particles[j].position]
    for particle in ourbox.particles[j].neighbourlist:
        nbpositionsj.append(particle.position)

    nbpositionsj = np.asarray(nbpositionsj)
    print(ourbox.particles[j].neighbourlist)
    ax.scatter(*np.asarray(ourbox.positions).T, color = "blue")
    ax.scatter(*nbpositionsj.T, color = "red")

    plt.show()









