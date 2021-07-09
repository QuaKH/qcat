#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import glob
import sys


def file_to_differential(path):
    s = ""
    with open(path, 'r') as reader:
        s = reader.read()[10:-3]
    
    rows = []
    cols = []
    vals = []
    
    bitmask = np.int64(2**32)
    
    for elem in s.split(", "):
        n = np.int64(elem)
        rows.append((abs(n) // bitmask) - 1)
        cols.append((abs(n) % bitmask) - 1)
        vals.append(1 if n > 0 else -1)
    
    return vals, (rows, cols)



def get_eigs(laplacian):
    # get smallest eigenpair

    eig_val_small = -1
    # if laplacian.det() != 0: # @TODO: somehow decide if this laplacian is nonsingular
    # eig_val_small, eig_vec_small = scipy.sparse.linalg.eigsh(laplacian, k=1, sigma=0, maxiter=1000)
    eig_vals = scipy.linalg.eigh(laplacian.toarray(), eigvals_only=True)
    
    # get only largest eigenvalue using sparse matrix
    # eig_val_big, eig_vec_big = scipy.sparse.linalg.eigsh(laplacian, k=1, maxiter=1000)

    # return (eig_val_small, eig_val_big)
    return eig_vals




def read_differential_files(dir):
    files = glob.glob(dir + "/*")

    dict = {}

    for file in files:
        data = file.split("_")
        i = int(data[-2])
        j = int(data[-1])
        
        dict[(i,j)] = file
    
    return dict



def get_diff_matrix(path):
    return scipy.sparse.coo_matrix(file_to_differential(path)).asfptype()

def write_laplacian_sparsity(laplacian, crossings, index, i, j):
    sparsity = laplacian.nnz

    with open("./laplacian_sparsity/knot_" + str(crossings) + "_" + str(index) + "_laspa", "a+") as writer:
        writer.write(str(i) + " " + str(j) + " " + str(laplacian.nnz) + " " + str(laplacian.shape[0]) + "\n")


def get_knot_eigs(dir, crossings, index):

    # @TODO: pass crossings and index to read_differential instead of parsing
    dict = read_differential_files(dir)
    keys = dict.keys()

    laplacians = {}

    for (i, j) in sorted(keys, reverse=True):

        d_i = get_diff_matrix(dict[(i,j)])

        laplacian = 0
        if (i-1,j) in keys:
            d_i_minus_1 = get_diff_matrix(dict[(i-1,j)])
            laplacian = (d_i_minus_1 * d_i_minus_1.transpose()) + (d_i.transpose() * d_i)
            laplacians[(i,j)] = laplacian
        else:
            laplacian = d_i.transpose() * d_i    
            laplacians[(i,j)] = laplacian

        # @TODO: save laplacian sparsity
        
        if (i+1,j) not in laplacians:
            laplacian = d_i * d_i.transpose()
            laplacians[i+1,j] = laplacian
    
    with open("./eigs/knot_" + str(crossings) + "_" + str(index) + "_eigs", "w+") as writer:
        for (i,j) in laplacians.keys():
            write_laplacian_sparsity(laplacians[(i,j)], crossings, index, i, j)

            if laplacians[(i,j)].shape == (1,1):
                continue

            # eig_small, eig_large = get_eigs(laplacians[(i,j)])
            eig_vals = get_eigs(laplacians[(i,j)])

            # write eigenvalues to file
            output_line = str(i) + " " + str(j) + " "
            writer.write(output_line)
            for e in eig_vals:
                writer.write(str(e) + " ")
            writer.write("\n")

if __name__ == "__main__":
    # get_knot_eigs(sys.argv[1], sys.argv[2], sys.argv[3])
    pass