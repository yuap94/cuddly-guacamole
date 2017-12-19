import numpy as np
#import system?


def enforce_pbc(r_vec, boxsize):
    for i, length in enumerate(boxsize):
        while r_vec[i] >= 0.5 * length:
            r_vec[i] -= length
        while r_vec[i] < -0.5 * length:
            r_vec[i] += length
    return r_vec

def verlet_neighbourlist(box, r_cut, r_skin):
	"""Verlet neighbourlist computation: for each of the particles in the array *particles*,
	we compute a list of particles that lie within a cutoff radius r_cut+r_skin. Returns 
	the array particles with an updated neighbourlist for each particle.

	arguments:
		box (Box object): includes array *box.particles* (array-like of Particle) listing all the
		Particle objects in the system
		r_cut: verlet cutoff radius
		r_skin: size of verlet skin region
	"""

	for i, particlei in enumerate(box.particles):
		particlei.neighbourlist = [] # clear neighbourlist for each particle		
	for particlej in box.particles[i+1:]:
		if np.linalg.norm(enforce_pbc(particlei.position - particlej.position,box.size)) < r_cut + r_skin
			particlei.neighbourlist.append(particlej) #this updates box.particles[i].neighbourlist as well?
			particlej.neighbourlist.append(particlei)

	return box

