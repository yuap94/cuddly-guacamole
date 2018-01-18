import numpy as np
import system
import numba as nb

@nb.jit
def verlet_neighbourlist(particles, r_cut, r_skin):
    """Verlet neighbourlist computation: for each of the particles in the array *box.particles*,
    we compute a neighbourlist of particles that lie within a cutoff radius r_cut+r_skin. 
    Returns an updated list of particles in the box with an updated neighbourlist for each particle.

    arguments:
        particles (list of Particle): listing of Particle objects in the system
        r_cut: verlet cutoff radius
        r_skin: size of verlet skin region
    """

    for i in range(len(particles)):
        particles[i].neighbourlist = [] # clear neighbourlist for each particle

    for i, particlei in enumerate(particles):
        for j, particlej in enumerate(particles[i+1:]):
            if np.linalg.norm(particlei.position - particlej.position) < r_cut + r_skin:
                particles[i].neighbourlist.append(particlej)
                particles[j+i+1].neighbourlist.append(particlei)

    return particles

