import numpy as np
import csv
import system
import input.inputgenerator
import metropolis
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


input.inputgenerator.gen_random_input(filename = "input/testinput.csv", n_particles= 3)


sysconfig = []

fid = open("input/testinput.csv","r")
fid_reader = csv.reader(fid, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
for row in fid_reader:
	sysconfig.append(row)

particles = []

# xenon: http://www.chegg.com/homework-help/questions-and-answers/2-values-lennard-jones-potential-xenon-found-177-kj-mol-410-respectively-values-lennard-jo-q14503605
# argon: http://pbx-brasil.com/outrasDisciplinas/DinMol/Notas/IIarea/aula203/papers/PhysRev.159.98.pdf
kb = 1.38064852*10**(-23) # boltzmann constant
sigma_argon = 3.405 # Å
epsilon_argon = 119.8 * kb * 10**(10) # N*Å

for i in range(len(sysconfig)):
	particles.append(system.Particle(position = np.asarray(sysconfig[i][0:3]), 
		charge = sysconfig[i][3], sigmaLJ = sigma_argon, epsilonLJ = epsilon_argon))

ourbox = system.Box(dimension = 3, size = np.array([1.0, 1.0, 1.0]), center = np.array([0.0, 0.0, 0.0]), 
					particles = particles, temp = 120.0)


width = sigma_argon / 10
n_steps = 10000
width = 1.0
n_skip = int(n_steps/100)
n_reuse_nblist = int(n_steps/50)
save_system_history = True
r_cut_LJ = 2.5 * 0.04 
r_skin_LJ = 2*n_reuse_nblist*width
ourbox, pos_history, pot_history = metropolis.mcmc(ourbox, n_steps, width, n_skip, n_reuse_nblist, save_system_history, r_cut_LJ, r_skin_LJ)

# print(ourbox.LJpotential)

# for i, pot in enumerate(pot_history):
# 	print(i)
# 	# print(pos)
# 	print(pot)


xyz = np.zeros((len(pos_history),3))
for i,pos_i in enumerate(pos_history):
 	xyz[i] = np.asarray(pos_i[0])
xyz = np.asarray(xyz)

fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')

ax.plot(*xyz.T, "-", color = "red") 
ax.scatter(*xyz.T, color = "blue")

ax = fig.add_subplot(122)

ax.plot(np.asarray(pot_history).T)

plt.show()


