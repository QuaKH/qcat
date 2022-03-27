#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import sqlite3
import glob
import threading
import queue
import pathlib


# In[3]:


# gets the file names where each differential matrix is stored (in base directory dir)
# returns dictionary with domain bigrading as keys and file name as values
def read_differential_files(dir):
    files = glob.glob(dir + "/*")

    d = {}

    for file in files:
        data = file.split("_")
        i = int(data[-2])
        j = int(data[-1])
        
        d[(i,j)] = file
    
    return d

# gets differential matrix data from file
def read_differential_from_file(path):
    vec = None
    shape = None
    with open(path, 'r') as reader:
        s = reader.readline().split(" ")
        shape = (int(s[0]), int(s[1]))
        vec = reader.readline()[10:-3]
    
    rows = []
    cols = []
    vals = []
    
    # bitmask to read rows and cols from Pari/GP sparse format
    bitmask = np.int64(2**32)
    
    for elem in vec.split(", "):
        n = np.int64(elem)
        rows.append((abs(n) // bitmask) - 1)
        cols.append((abs(n) % bitmask) - 1)
        vals.append(1 if n > 0 else -1)
    
    return (vals, (rows, cols)), shape

# returns differential matrix stored at file path as scipy sparse matrix
def get_diff_matrix(path):
    (vals, (rows, cols)), s = read_differential_from_file(path)
    return scipy.sparse.coo_matrix((vals, (rows, cols)), shape=s).asfptype()


# In[4]:


# computes smallest and largest eigenvalues of laplacian matrix
def get_lap_eigs(laplacian, num_zero_eigs):
    # handle special case where laplacian is 1x1
    if laplacian.shape[0] == 1:
        val = laplacian.toarray()[0][0]
        return val, val


    largest = scipy.sparse.linalg.eigsh(laplacian, which="LM", return_eigenvectors=False, k=1, tol=10e-5)[0]
    
    # Only one nonzero eigenvalue
    if num_zero_eigs == laplacian.shape[0] - 1:
        return largest, largest
    
    smallest = scipy.sparse.linalg.eigsh(laplacian, which="SM", return_eigenvectors=False, k=num_zero_eigs + 1, tol=10e-5)[0]

    return smallest, largest


# In[5]:


# compute the dimension of the homology indicates the number of eigenvalues of laplacian
# that are exactly zero
# pass None for differential if zero matrix
def get_num_zero_eigs(d_i_minus_1, d_i):
    img_dim = 0
    ker_dim = 0
    if (d_i != None):
        img_dim = np.linalg.matrix_rank(d_i.toarray())
        ker_dim = d_i.shape[1]

    if (d_i_minus_1 != None):
        ker_dim = d_i_minus_1.shape[0] - np.linalg.matrix_rank(d_i_minus_1.toarray())
    
    return ker_dim - img_dim


# In[6]:


# converts pd code string from format:
# 1, 2, 3, 4; 1, 2, 3, 4; ...
# to
# [[1,2,3,4],[1,2,3,4]...]
def format_pd_code(pd_code_string):
    columns = pd_code_string.split(";")
    columns = map(lambda s: s.strip(), columns)
    return "[ [" + "], [".join(columns) + "] ]"


# In[7]:


# creates a dictionary of laplacians accesible by (i,j) grading
def get_laplacian_dict(differential_dict):
    keys = differential_dict.keys()
    laplacians = {}

    for (i, j) in sorted(keys, reverse=True):
        d_i = get_diff_matrix(differential_dict[(i,j)])

        if (i-1,j) in keys:
            d_i_minus_1 = get_diff_matrix(differential_dict[(i-1,j)])
            laplacian = (d_i_minus_1 * d_i_minus_1.transpose()) + (d_i.transpose() * d_i)         
            laplacians[(i,j)] = (laplacian, get_num_zero_eigs(d_i_minus_1, d_i))
        else:
            laplacian = d_i.transpose() * d_i    
            laplacians[(i,j)] = (laplacian, get_num_zero_eigs(None, d_i))
        
        if (i+1,j) not in laplacians:
            laplacian = d_i * d_i.transpose()
            laplacians[i+1,j] = (laplacian, get_num_zero_eigs(d_i, None))
    
    return laplacians

# computes eigenvalues of all laplacians for specified knot type and writes to file
# dir specifies directory with differential matrices for this knot diagram
def get_knot_eigs(dir, pd_code, crossings, db_path, db_table_name):
    differential_dict = read_differential_files(dir)
    laplacians = get_laplacian_dict(differential_dict)

    con = sqlite3.connect(db_path)
    cur = con.cursor()


    for (i,j) in laplacians.keys():

        laplacian = laplacians[(i,j)][0]
        num_zero_eigs = laplacians[(i,j)][1]

        eig_vals = get_lap_eigs(laplacian, num_zero_eigs)

        # (pd-code, crossings, i, j, num-zero-eigs, laplacian-dim, laplacian-nnz, min-eig, max-eig)
        entry = (format_pd_code(pd_code), crossings, i, j, np.shape(laplacian)[0], laplacian.nnz, int(num_zero_eigs), eig_vals[0], eig_vals[1])
        
        cur.execute("INSERT INTO " + db_table_name + " VALUES (?,?,?,?,?,?,?,?,?)", entry)

    con.commit()
    con.close()

if __name__ == "__main__":
    import sys
    get_knot_eigs(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5])
