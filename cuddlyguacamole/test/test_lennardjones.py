import numpy as np
import numpy.testing as npt
import lennardjones
import system

def test_LJ_potential_ij(dim):


	kb = 1.38064852*10**(-23) # boltzmann constant

	sigma_argon = 3.405 # Å
	epsilon_argon = 119.8 * kb * 10**(10) # N*Å

	r = np.linalg.norm(np.random.rand(dim))
	r_c = 2.5*sigma_argon
	n_skip = 10
	width = r_c / 10
	r_skin = 2*n_skip*width

	r_cut = r_c + r_skin

	LJpot_ArAr_ref = (4*epsilon_argon*((sigma_argon/r)**(12)-(sigma_argon/r)**(6)) 
					- 4*epsilon_argon*((sigma_argon/r_cut)**(12)-(sigma_argon/r_cut)**(6)))

	LJpot_ArAr = lennardjones.LJ_potential_ij(r, sigma_argon, epsilon_argon, sigma_argon, epsilon_argon, r_c, r_skin)
	npt.assert_almost_equal(LJpot_ArAr, LJpot_ArAr_ref)



def test_LJ_potential_ij(dim):

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

	r_Ar1Ar2 = np.linalg.norm(argon_1.position - argon_2.position)
	r_Ar1Ar3 = np.linalg.norm(argon_1.position - argon_3.position)
	r_Ar1Xe1 = np.linalg.norm(argon_1.position - xenon_1.position)
	r_Ar1Xe2 = np.linalg.norm(argon_1.position - xenon_2.position)
	r_Ar1Xe3 = np.linalg.norm(argon_1.position - xenon_3.position)

	r_Ar2Ar3 = np.linalg.norm(argon_2.position - argon_3.position)
	r_Ar2Xe1 = np.linalg.norm(argon_2.position - xenon_1.position)
	r_Ar2Xe2 = np.linalg.norm(argon_2.position - xenon_2.position)
	r_Ar2Xe3 = np.linalg.norm(argon_2.position - xenon_3.position)
	
	r_Ar3Xe1 = np.linalg.norm(argon_3.position - xenon_1.position)
	r_Ar3Xe2 = np.linalg.norm(argon_3.position - xenon_2.position)
	r_Ar3Xe3 = np.linalg.norm(argon_3.position - xenon_3.position)

	r_Xe1Xe2 = np.linalg.norm(xenon_1.position - xenon_2.position)
	r_Xe1Xe3 = np.linalg.norm(xenon_1.position - xenon_3.position)
	
	r_Xe2Xe3 = np.linalg.norm(xenon_2.position - xenon_3.position)

	# neighbourlists.....
	LJ_pot_ref = (lennardjones.LJ_potential_ij(r_Ar1Ar2, sigma_argon, epsilon_argon, sigma_argon, epsilon_argon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Ar1Ar3, sigma_argon, epsilon_argon, sigma_argon, epsilon_argon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Ar1Xe1, sigma_argon, epsilon_argon, sigma_xenon, epsilon_xenon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Ar1Xe2, sigma_argon, epsilon_argon, sigma_xenon, epsilon_xenon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Ar1Xe3, sigma_argon, epsilon_argon, sigma_xenon, epsilon_xenon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Ar2Ar3, sigma_argon, epsilon_argon, sigma_argon, epsilon_argon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Ar2Xe1, sigma_argon, epsilon_argon, sigma_xenon, epsilon_xenon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Ar2Xe2, sigma_argon, epsilon_argon, sigma_xenon, epsilon_xenon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Ar2Xe3, sigma_argon, epsilon_argon, sigma_xenon, epsilon_xenon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Ar3Xe1, sigma_argon, epsilon_argon, sigma_xenon, epsilon_xenon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Ar3Xe2, sigma_argon, epsilon_argon, sigma_xenon, epsilon_xenon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Ar3Xe3, sigma_argon, epsilon_argon, sigma_xenon, epsilon_xenon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Xe1Xe2, sigma_xenon, epsilon_xenon, sigma_xenon, epsilon_xenon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Xe1Xe3, sigma_xenon, epsilon_xenon, sigma_xenon, epsilon_xenon, r_c, r_s)
				 +lennardjones.LJ_potential_ij(r_Xe2Xe3, sigma_xenon, epsilon_xenon, sigma_xenon, epsilon_xenon, r_c, r_s))

	ourbox.compute_LJ_potential(r_c, r_s)
	npt.assert_almost_equal(LJ_pot_ref, ourbox.LJpotential)



