import numpy as np
import numpy.testing as npt
import neighbourlist
import system

def test_verlet_neighbourlist(dim):

	kb = 1.38064852*10**(-23) # boltzmann constant

	sigma_argon = 3.405 # Å
	epsilon_argon = 119.8 * kb * 10**(10) # N*Å

	sigma_xenon = 4.07 # Å
	epsilon_xenon = 225.3 * kb * 10**(10) # Å

	r_c = 2.5*0.5*(sigma_argon+sigma_xenon)
	width = r_c / 10
	n_skip = 10
	r_s = 2*n_skip*width

	boxsize = np.array([1.0, 1.0, 1.0])

	pos1 = system.enforce_pbc(np.random.rand(dim), boxsize)
	pos2 = system.enforce_pbc(np.random.rand(dim), boxsize)
	pos3 = system.enforce_pbc(np.random.rand(dim), boxsize)
	pos4 = system.enforce_pbc(np.random.rand(dim), boxsize)
	pos5 = system.enforce_pbc(np.random.rand(dim), boxsize)
	pos6 = system.enforce_pbc(np.random.rand(dim), boxsize)

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
		if not particles[i].neighbourlist == ourbox.particles[i].neighbourlist:
			print("neighbourlists incorrect")
		# print("particles[" + repr(i) + "] = " + repr(particles[i].neighbourlist))
		# print("ourbox.particles[" + repr(i) + "] = " + repr(ourbox.particles[i].neighbourlist))











