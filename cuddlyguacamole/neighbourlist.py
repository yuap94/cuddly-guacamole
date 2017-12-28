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
	"""Verlet neighbourlist computation: for each of the particles in the array *box.particles*,
	we compute a neighbourlist of particles that lie within a cutoff radius r_cut+r_skin. 
	Returns the box with an updated neighbourlist for each particle.

	arguments:
		box (Box object): includes array *box.particles* (array-like of Particle) listing all the
		Particle objects in the system
		r_cut: verlet cutoff radius
		r_skin: size of verlet skin region
	"""

	for i in range(len(box.particles)):
		box.particles[i].neighbourlist = [] # clear neighbourlist for each particle

	for i, particlei in enumerate(box.particles):		
		for particlej in box.particles[i+1:]:
			if np.linalg.norm(enforce_pbc(particlei.position - particlej.position, box.size)) < r_cut + r_skin
				box.particles[i].neighbourlist.append(box.particles[j])
				box.particles[j].neighbourlist.append(box.particles[i])

	return box

