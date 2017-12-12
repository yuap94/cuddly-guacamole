import numpy as np


############################################################################
# Define Particle and Box classes 
############################################################################
class Particle(object):
    """A particle. Particles have the following properties:

    Attributes:
        position: a one to three dimensional position vector of numpy.array type
        charge (float): charge of the particle
    """

    def __init__(self, position, charge, mass = 0.0):
        """Return a Particle object whose position is *position*,
        electric charge is *charge*, and mass is *mass* (can be left out
        if mass is not relevant, default is 0.0)"""
        self.position = position
        self.charge = charge
        self.mass = mass


class Box(object):
	"""A 1 to 3-dimensional square box in space, containing charged particles:

	Attributes:
		dimension (int): dimension of the box (1-3)
	    size (float): size of the (quadratic) box
	    center (float): 1-3d numpy array locating the center of the box
	    particles (Particle): numpy array of particles in the box (particle objects 
	    have charge, position and possibly mass)
	"""

	def __init__(self, dimension, size, center, particles):
	    """Return a Box object of dimension *dimension* (between 1 and 3),
	    whose length(&breadth&height) is *size*, is centered at *center*, 
	    and contains the particles in the numpy array *particles*"""
	    self.dimension = dimension
	    self.size = size
	    self.center = center
	    self.particles = particles

#############################################################################