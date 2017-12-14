import numpy as np

def verlet_neighbourlist(particles, r_cut, r_skin):
	
	"""Verlet neighbourlist computation: for each of the particles in the array *particles*,
	we compute a list of particles (among those in the array *particles*) that lie within
	a cutoff radius r_cut+r_skin. Returns the array particles with an updated neighbourlist for
	each particle.
    
    arguments:
    	particles (array-like of Particle): list of all the Particle objects in the system
    	r_cut: verlet cutoff radius (mins r_skin)
    	r_skin: size of skin region
    """

	for particle in particles:
		particle.neighbourlist = [] # clear neighbourlist for all particles

	for i, particlei in enumerate(particles):
		for particlej in particles[i+1:]:
			if np.linalg.norm(particlei.position - particlej.position) < r_cut + r_skin
				particlei.neighbourlist.append(particlej)
				particlej.neighbourlist.append(particlei)

