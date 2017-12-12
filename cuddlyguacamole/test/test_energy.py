import numpy as np
import numpy.testing as npt
from ..energy import pme
from ..system import Particle
from ..system import Box

	def test_LR_energy():
		a = np.array([Particle(np.array([1, 2, 3.0]), -1.0), Particle(np.array([1, 2, 3.0]), 1.0)]) 
		#...
		box = Box(a,...)

		system = np.array(box...)

		npt.assert_allclose(pme(system, 0.1), energycomputedbyhand, atol=0.001)


