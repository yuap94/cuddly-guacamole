import test.test_lennardjones
import test.test_neighbourlist
import test.test_metropolis

test.test_lennardjones.test_LJ_potential_ij(3)

test.test_lennardjones.test_LJ_potential(3)

test.test_neighbourlist.test_verlet_neighbourlist(3)

# test_neighbourlist.plot_neighbourlist()

# n_same, n_diff = 0, 0
# for i in range(100):
#     ourbox_ref, ourbox = test_metropolis.test_mcmc_step(3, True)
#     if (ourbox_ref.positions[0] == ourbox.positions[0]).all():
#         n_same += 1
#     else:
#         n_diff += 1
# print("n_same = " + repr(n_same) + ", n_diff = " + repr(n_diff))

n_steps = 1000
n_reuse_nblist = int(n_steps/100)
test.test_metropolis.test_mcmc(3, n_steps, n_skip = int(n_steps/10), n_reuse_nblist = n_reuse_nblist)