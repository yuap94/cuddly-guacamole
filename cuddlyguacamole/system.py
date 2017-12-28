import numpy as np
import lennardjones
import neighbourlist


############################################################################
# Define Particle and Box classes 
############################################################################
class Particle(object):
    """A particle. Particles have the following properties:

    Attributes:
        position: a one to three dimensional position vector of numpy.array type
        charge (float): charge of the particle
        self.sigmaLJ (float): distance at which LJ potential is 0 for particle type
        self.epsilonLJ (float): depth of LJ potential well for particle type
        neighbourlist (array like of Particle): a list of of neighbours of the particle
        (verlet neighbourlist of particles within cutoff+skin radius) 
    """

    def __init__(self, position, charge, sigmaLJ = 1.0, epsilonLJ = 1.0):
        """Return a Particle object whose position is *position*,
        electric charge is *charge*, sigma and epsilon for the Leonard-Jones
        potential is *sigmaLJ* and *epsilonLJ*, and neighbourlist is *neighbourlist*"""
        self.position = position
        self.charge = charge
        self.sigmaLJ = sigmaLJ
        self.epsilonLJ = epsilonLJ
        self.neighbourlist = None


class Box(object):
	"""A 1 to 3-dimensional square box in space, containing charged particles:

	Attributes:
		dimension (int): dimension of the box (1-3)
	    size (*dimension*-dimensional numpy array of float): 1d-3d numpy array giving size of the box in each direction
	    center (*dimension*-dimensional numpy array of loat): 1d-3d numpy array locating the center of the box
	    particles (*dimension*-dimensional numpy array of Particle): 1d-3d numpy array of particles in the box (particle objects 
	    		have charge, position and possibly mass)
	    positions (python list of *dimension*-dimensional numpy arrays with the position of all the particles, for ease of use... only
	    		create if needed using a set method? and always store if created?)
	    		region (*dimension*x2-dimensional numpy array of float): specifies the region
	    		in space that the box covers (automatically computed from the center and the size of the box?)
	    LJpotential (float): Lennard Jones Potential of the system (calculated based on the positions of the particles in *particles*)
	    temp (float): temperature in the box
	    LJneighbourlist (list of numpy arrays of int): a list with of same size as *particles*,
	    											 with a numpy array of int for each particle listing the indices of the 
	    											 neighbouring particles (particles within the LJ cutoff radius)
		r_c_LJ (float): cutoff radius for LJ potential calculation
		r_skin_LJ (float): size of skin region for LJ potential calculation 
	"""

	def __init__(self, dimension, size, center, particles, temp=273.15, LJneighbourlist = []):
	    """Return a Box object of dimension *dimension* (between 1 and 3),
	    whose length(&breadth&height) is *size*, is centered at *center*, 
	    and contains the particles in the numpy array *particles*"""
	    self.dimension = dimension
	    self.size = size
	    self.center = center
	    self.particles = particles
	    self.region = np.transpose(np.array([center-size/2, center + size/2]))
	    self.LJpotential = None
	    self.positions = None
	    self.temp = temp
	    self.LJneighbourlist = LJneighbourlist

	@LJneighbourlist.setter    # why use decorator??
    def compute_LJneighbourlist(self, r_cut, r_skin):
    	self.LJneighbourlist = neighbourlist.LJneighbourlist(self, r_cut, r_skin)

	@LJpotential.setter    # why use decorator??
    def compute_LJ_potential(self, r_cut, r_skin):
    	self.LJpotential = lennardjones.LJ_potential(self, r_cut, r_skin)

	@positions.setter # why use decorator??
	def make_positions_list(self):
		self.positions = [] # store initial position for each particle in list
    	for particle in self.particles:
        	self.positions.append(particle.position)



#############################################################################