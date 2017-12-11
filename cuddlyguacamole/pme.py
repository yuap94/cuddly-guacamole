
%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
from time import time


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

    # def withdraw(self, amount):
    #     """Return the balance remaining after withdrawing *amount*
    #     dollars."""
    #     if amount > self.balance:
    #         raise RuntimeError('Amount greater than available balance.')
    #     self.balance -= amount
    #     return self.balance


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

class Ensemble(object):
	"""A 1 to 3-dimensional space, containing boxes of charged particles:

	Attributes:
		dimension (int): dimension of the space
	    extension (float): total size of the (quadratic) space
	    boxes (numpy array of Box): array of Box objects in the space
	"""

	def __init__(self, dimension, size, box, box_position, no_of_boxes):
	    """Return an ensemble object whose length(&breadth&height) is *size*,
	    is centered at *center*, and is divided into the Box objects in the numpy
	    array *boxes*"""
	    self.dimension = dimension	  
	    self.size = size
	    self.box = box
	    self.no_of_boxes = no_of_boxes
	    self.box_position = box_position
	    self.boxes = np.array[(.....)] # construct the array of boxes in the ensemble using the information about
	    # the box size and the amount of boxes/position of boxes


def spreading_fcn(d, r, sigma):
    """Return value of d-dimensional gaussian spreading function with std sigma at point r"""
    return np.exp(-np.inner(r,r)/(2*sigma**2))/(2*np.pi*sigma**2)**(3/2)



def PME(ensemble, sigma): 
    '''Particle mesh Ewald summation, computes the long range
    Ewald (coloumb) interaction energy of an ensemble of charged particles.
    
    arguments:
    	ensemble (object of class Ensemble): attributes of the ensemble includes an its extension in space
    	an array of box objects, each with attributes extension (giving their size) and center (giving the 
    	position of the center of the box).

    	The boxes in turn have a numpy array of particles, (i.e. objects of the class Particle).
        box: a box object, giving the dimensions of the box the particles are in
        space: a space object, specifying the dimensions of the total space we're considering,
        how many boxes are in our space and where they are located (their centers and dimensions)
    '''
	
	# First we compute the LR charge density at the center of each box, using the charges and positions
	# of the particles in the box
	M = ensemble.boxes.size
	rho_L = np.zeros(M) #
	d = ensemble.dimension
	L = box.size

	for m in range(0,M):
		for box in ensemble.boxes:
			n = box.center
			for j in range(size(box.particles)):
				rho_L[m] += box.particles[j].charge*spreading_fcn(d, particles[m].position - particles[j].position + n*L, sigma)  


	# r = np.zeros((steps,2))
    # r[0] = np.random.randn(2)

    # for i in range(steps-1):
        
    #     r[i+1] = r[i] + dt/m*p[i] - dt**2/(2*m)*pot_grad(r[i])

    return r, p