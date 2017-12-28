import numpy as np
import csv
import system
import input.inputgenerator


input.inputgenerator.gen_random_input("input/testinput.csv", 3)


sysconfig = []

fid = open("input/testinput.csv","r")
fid_reader = csv.reader(fid, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
for row in fid_reader:
	sysconfig.append(row)

particles = []

for i in range(len(sysconfig)):
	particles.append(system.Particle(position = sysconfig[i][0:3], 
		charge = sysconfig[i][3], sigmaLJ = 1.0, epsilonLJ = 1.0))

ourbox = system.Box(dimension = 3, size = np.array([1.0, 1.0, 1.0]), center = np.array([0.0, 0.0, 0.0]), 
					particles = particles, temp = 200.0)

print(ourbox)