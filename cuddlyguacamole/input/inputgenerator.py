import csv
import numpy as np
import pbc


def gen_random_input_3D(filename, n_particles, boxsize, r_c):

    fid = open(filename, "w")
    fid_writer = csv.writer(fid, delimiter=',')

    sysconfig = np.zeros((n_particles,4)) # generate 3d particle position vectors and particle charge
    sysconfig.T[3] = np.random.randn(1, n_particles) # the charge
    sysconfig[0][0:3] = pbc.enforce_pbc(np.random.randn(3) * boxsize / 4, boxsize) # position of particle0
    for i in range(1,len(sysconfig)):
        sysconfig[i][0:3] = pbc.enforce_pbc(sysconfig[i-1][0:3] + np.random.randn(3) * r_c / (n_particles*8), boxsize) # perturb position of particle0 to get initial pos of particle1
    sysconfig = sysconfig.tolist()

    # print(sysconfig)

    for i in range(len(sysconfig)):
        fid_writer.writerow(sysconfig[i])

    fid.close()

# def gen_random_input_3D(filename, n_particles):

#     fid = open(filename, "w")
#     fid_writer = csv.writer(fid, delimiter=',')

#     sysconfig = (np.random.randn(n_particles,4))/np.array([10, 10, 10, 1]) # divide by 10 to keep particles close together early on
#     sysconfig = sysconfig.tolist()

#     for i in range(len(sysconfig)):
#         fid_writer.writerow(sysconfig[i])

#     fid.close()
