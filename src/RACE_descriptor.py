import os
import re
import numpy as np
import pandas as pd
from ase.io import read
from ase.neighborlist import NeighborList
from ase.data import atomic_numbers

STRUCTURE_FOLDER = "structures"
OUTPUT_FILE = "descriptor_RACE.csv"

# =========================
# parâmetros do descritor
# =========================

R_CUT = 4.5

N_R = 16
SIGMA_R = 0.08
R_MIN = 1.0

N_A = 16
SIGMA_A_DEG = 7.0

NORMALIZE_LOCAL = True
NORMALIZE_GLOBAL = True

EPS = 1e-12


# =========================
# utilidades
# =========================

def numeric_sort_key(filename):
    m = re.match(r"(\d+)", filename)
    return int(m.group(1)) if m else float("inf")


def smooth_cutoff(r):

    if r >= R_CUT:
        return 0.0

    return 0.5 * (np.cos(np.pi * r / R_CUT) + 1.0)


def chem_weight(Zi, Zj):

    return np.sqrt(Zi * Zj)


def safe_l1_normalize(x):

    s = np.sum(x)

    if s > EPS:
        return x / s

    return x


# =========================
# ambiente local
# =========================

def compute_local_environment(atoms, i, Z, nl, r_centers, a_centers, sigma_a):

    Zi = Z[i]

    pos = atoms.get_positions()
    cell = atoms.get_cell()

    rdf_i = np.zeros(N_R)
    adf_i = np.zeros(N_A)
    F_i = np.zeros((N_R, N_A))

    idxs, offs = nl.get_neighbors(i)

    neigh_vecs = []
    neigh_r = []
    neigh_w = []

    for j, off in zip(idxs, offs):

        rj = atoms.positions[j] + np.dot(off, cell)

        rij_vec = rj - pos[i]
        rij = np.linalg.norm(rij_vec)

        if rij < EPS:
            continue

        wc = smooth_cutoff(rij)

        if wc <= 0:
            continue

        Zj = Z[j]

        w = wc * chem_weight(Zi, Zj)

        neigh_vecs.append(rij_vec)
        neigh_r.append(rij)
        neigh_w.append(w)

        rdf_i += w * np.exp(-0.5 * ((rij - r_centers) / SIGMA_R) ** 2)

    m = len(neigh_vecs)

    if m >= 2:

        neigh_vecs = np.array(neigh_vecs)
        neigh_r = np.array(neigh_r)
        neigh_w = np.array(neigh_w)

        Kr = np.exp(-0.5 * ((neigh_r[:, None] - r_centers) / SIGMA_R) ** 2)

        for p in range(m):

            v1 = neigh_vecs[p]
            n1 = np.linalg.norm(v1)

            if n1 < EPS:
                continue

            for q in range(p + 1, m):

                v2 = neigh_vecs[q]
                n2 = np.linalg.norm(v2)

                if n2 < EPS:
                    continue

                cosang = np.dot(v1, v2) / (n1 * n2)
                cosang = np.clip(cosang, -1, 1)

                theta = np.arccos(cosang)

                wa = np.exp(-0.5 * ((theta - a_centers) / sigma_a) ** 2)

                w_pair = neigh_w[p] * neigh_w[q]

                adf_i += w_pair * wa

                radial_sum = Kr[p] + Kr[q]

                F_i += (w_pair * radial_sum)[:, None] * wa[None, :]

    if NORMALIZE_LOCAL:

        rdf_i = safe_l1_normalize(rdf_i)
        adf_i = safe_l1_normalize(adf_i)
        F_i = safe_l1_normalize(F_i)

    return rdf_i, adf_i, F_i


# =========================
# agregação
# =========================

def aggregate_descriptors(rdf_list, adf_list, F_list):

    rdf_stack = np.array(rdf_list)
    adf_stack = np.array(adf_list)
    F_stack = np.array(F_list)

    rdf_mean = rdf_stack.mean(axis=0)
    rdf_std = rdf_stack.std(axis=0)

    adf_mean = adf_stack.mean(axis=0)
    adf_std = adf_stack.std(axis=0)

    F_mean = F_stack.mean(axis=0)
    F_std = F_stack.std(axis=0)

    if NORMALIZE_GLOBAL:

        rdf_mean = safe_l1_normalize(rdf_mean)
        adf_mean = safe_l1_normalize(adf_mean)
        F_mean = safe_l1_normalize(F_mean)

    return rdf_mean, rdf_std, adf_mean, adf_std, F_mean, F_std


# =========================
# descritor principal
# =========================

def compute_descriptor(atoms):

    symbols = atoms.get_chemical_symbols()
    Z = np.array([atomic_numbers[s] for s in symbols])

    n = len(atoms)

    cutoffs = [R_CUT / 2] * n

    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    r_centers = np.linspace(R_MIN, R_CUT, N_R)

    a_centers = np.deg2rad(np.linspace(0, 180, N_A))
    sigma_a = np.deg2rad(SIGMA_A_DEG)

    rdf_list = []
    adf_list = []
    F_list = []

    for i in range(n):

        rdf_i, adf_i, F_i = compute_local_environment(
            atoms, i, Z, nl, r_centers, a_centers, sigma_a
        )

        rdf_list.append(rdf_i)
        adf_list.append(adf_i)
        F_list.append(F_i)

    rdf_mean, rdf_std, adf_mean, adf_std, F_mean, F_std = aggregate_descriptors(
        rdf_list, adf_list, F_list
    )

    feats = {}

    for k in range(N_R):

        feats[f"rdf_mean_{k:02d}"] = rdf_mean[k]
        feats[f"rdf_std_{k:02d}"] = rdf_std[k]

    for k in range(N_A):

        feats[f"adf_mean_{k:02d}"] = adf_mean[k]
        feats[f"adf_std_{k:02d}"] = adf_std[k]

    flat_mean = F_mean.reshape(-1)
    flat_std = F_std.reshape(-1)

    for t in range(N_R * N_A):

        feats[f"joint_mean_{t:04d}"] = flat_mean[t]
        feats[f"joint_std_{t:04d}"] = flat_std[t]

    # estatísticas químicas globais

    feats["n_atoms"] = n
    feats["z_mean"] = np.mean(Z)
    feats["z_std"] = np.std(Z)
    feats["z_min"] = np.min(Z)
    feats["z_max"] = np.max(Z)

    return feats


# =========================
# pipeline
# =========================

def main():

    files = [f for f in os.listdir(STRUCTURE_FOLDER) if f.endswith(".cif")]

    files.sort(key=numeric_sort_key)

    rows = []

    for f in files:

        print("Processing", f)

        atoms = read(os.path.join(STRUCTURE_FOLDER, f))

        feats = compute_descriptor(atoms)

        feats["filename"] = f

        rows.append(feats)

    df = pd.DataFrame(rows)

    cols = ["filename"] + [c for c in df.columns if c != "filename"]

    df = df[cols]

    df.to_csv(OUTPUT_FILE, index=False)

    print("Saved", OUTPUT_FILE)


if __name__ == "__main__":

    main()
