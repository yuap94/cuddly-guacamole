import csv
import numpy as np


def gen_random_input(filename, n_particles):

	fid = open(filename, "w")
	fid_writer = csv.writer(fid,delimiter=',')

	sysconfig = (np.random.rand(n_particles,4) - np.array([0.5,0.5,0.5,0.5])) * np.array([1,1,1,2])
	sysconfig = sysconfig.tolist()

	for i in range(len(sysconfig)):
		fid_writer.writerow(sysconfig[i])

	fid.close()
