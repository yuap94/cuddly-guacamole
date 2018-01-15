def enforce_pbc(r_vec, boxsize):
    for i, length in enumerate(boxsize):
        while r_vec[i] >= 0.5 * length:
            r_vec[i] -= length
        while r_vec[i] < -0.5 * length:
            r_vec[i] += length
    return r_vec