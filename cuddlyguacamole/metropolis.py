import numpy as np
import lennardjones
import energy

def mcmc_step(r_vec, box, e_vec=None, width=0.2):
    if e_vec is None:
        e_vec = potential(r_vec, box)
    r_trial = r_vec + width * (np.random.rand(*r_vec.shape) - 0.5)
    e_trial = potential(r_trial, box)
    if e_trial < e_vec or np.random.rand() < np.exp(e_vec - e_trial):
        return r_trial, e_trial
    return r_vec, e_vec

def mcmc(r_vec_init, n_steps, box, width=0.2, n_skip=1):
    r_vec = [np.asarray(r_vec_init)]
    e_vec = [potential(r_vec[-1], box)]
    r, e = r_vec[0], e_vec[0]
    for i in range(n_steps):
        for j in range(n_skip):
            r, e = mcmc_step(r, box, e_vec=e, width=width)
        r_vec.append(r)
        e_vec.append(e)
    return np.asarray(r_vec), np.asarray(e_vec)

r, e = mcmc(r_opt[-1], 10000, [6.0, 6.0], n_skip=100)

fig, ax = plt.subplots(figsize=(6, 6))
for i in range(3):
    ax.plot(*r[:, i, :].T, '--', color='grey', alpha=0.1)
    ax.scatter(*r[:, i, :].T, c=np.arange(r.shape[0]))
ax.set_aspect('equal')