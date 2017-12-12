
%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
from time import time
from system import Particle, Box


##### Compute LR energy with PME:

def spreading_fcn(d, r, sigma):
    """Return value of d-dimensional gaussian spreading function with std sigma at point r"""
    return np.exp(-np.inner(r,r)/(2*sigma**2))/(2*np.pi*sigma**2)**(3/2)


def interpolate(phi_L_discrete,...):
	...
	return phi_L_continuous

def PME(system, sigma): 
    '''Particle mesh Ewald summation, computes the long range
    Ewald (coloumb) interaction energy of an ensemble of charged particles.
    
    arguments:
    	system (array like of Box?): array of boxes  (each with a center, size and array particles of type Particle)
    	sigma: std of the gaussian spreading functions...
    '''
	
	# First we compute the LR charge density at the center of each box, using the charges and positions
	# of the particles in the box
	M = system.boxes[1].size
	rho_L = np.zeros(M) #
	d = system.dimension
	L = box.size

	for m in range(0,M):
		for box in system.boxes:
			n = box.center
			for j in range(size(box.particles)):
				rho_L[m] += box.particles[j].charge*spreading_fcn(d, particles[m].position - particles[j].position + n*L, sigma)  


	rho_L_hat = FFT(rho_L)

	for k in range(rho_hat.size):
		phi_L_hat[k] = rho_L/(epsilon0*k^2)

	phi_L[m] = IFFT(phi_L_hat)

	system.LRpotential = interpolate(phi_L)

    return system