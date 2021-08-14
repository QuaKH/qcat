#!/usr/bin/env python
# coding: utf-8

# In[67]:


import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import glob
import threading
import queue
import pathlib
import psutil
import torch

# In[61]:


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


# In[4]:


def get_lap_eigs(laplacian, num_zero_eigs):
    smallest, largest = None, None
    if torch.cuda.is_available():
        device = torch.device("cuda")

        i = torch.LongTensor(np.vstack((laplacian.row, laplacian.col)))
        v = torch.FloatTensor(laplacian.data)
        matrix = torch.sparse.FloatTensor(i, v, torch.Size(laplacian.shape))

        eigs = torch.linalg.eigvals(matrix)
        eigs_list = np.sort(eigs.numpy())
        return eigs_list[num_zero_eigs], eigs_list[-1]
    else:
        smallest = scipy.sparse.linalg.eigsh(laplacian, which="SM", return_eigenvalues=False, k=num_zero_eigs + 1, tol=10e-5)[0]
        largest = scipy.sparse.linalg.eigsh(laplacian, which="LM", return_eigenvalues=False, k=1, tol=10e-5)[0]
    return smallest, largest


# In[5]:


def read_differential_files(dir):
    files = glob.glob(dir + "/*")

    d = {}

    for file in files:
        data = file.split("_")
        i = int(data[-2])
        j = int(data[-1])
        
        d[(i,j)] = file
    
    return d


# In[56]:


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

# get_num_zero_eigs(get_diff_matrix("./KhoHo/differentials/knot_6_1/knot_6_1__d_-1_1"), get_diff_matrix("./KhoHo/differentials/knot_6_1/knot_6_1__d_0_1"))


# In[62]:



def get_diff_matrix(path):
    (vals, (rows, cols)), s = read_differential_from_file(path)
    return scipy.sparse.coo_matrix((vals, (rows, cols)), shape=s).asfptype()

def write_laplacian_sparsity(laplacian, crossings, index, i, j):
    sparsity = laplacian.nnz

    with open("./laplacian_sparsity/knot_" + str(crossings) + "_" + str(index) + "_laspa", "a+") as writer:
        writer.write(str(i) + " " + str(j) + " " + str(laplacian.nnz) + " " + str(laplacian.shape[0]) + "\n")

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

def get_knot_eigs(dir, crossings, index):
    differential_dict = read_differential_files(dir)
    laplacians = get_laplacian_dict(differential_dict)
    
    with open("./eigs/knot_" + str(crossings) + "_" + str(index) + "_eigs", "w+") as writer:
        for (i,j) in laplacians.keys():
            # record sparsity of laplacian
            write_laplacian_sparsity(laplacians[(i,j)][0], crossings, index, i, j)

            # get eigenvalues of laplacian
            smallest, largest = get_lap_eigs(laplacians[(i,j)][0], laplacians[(i,j)][1])

            # write eigenvalues to file
            output_line = str(i) + " " + str(j) + " " + str(smallest) + " " + str(largest)
            writer.write(output_line)
            writer.write("\n")

# get_knot_eigs("./KhoHo/differentials/knot_7_1/", 7, 1)


# In[ ]:


q = queue.Queue()

def worker():
    while (True):
        task = q.get()
        # print("working on " + str(task["crossings"]) + "_" + str(task["index"]))
        get_knot_eigs(task["path"], task["crossings"], task["index"])
        
        q.task_done()

def run_prague():
    MAX_THREAD_NUM = psutil.cpu_count()
    threads = [threading.Thread(target=worker, daemon=True).start() for i in range(MAX_THREAD_NUM)]

    index_count = [1, 1, 2, 3, 7, 21, 49, 165]
    for crossings in range(3,11):
        for index in range(1,index_count[crossings - 3] + 1):
            path = "./KhoHo/differentials/knot_" + str(crossings) + "_" + str(index)
            dir = pathlib.Path(path)
            
            if dir.is_dir():
                q.put({"path":path, "crossings": crossings, "index": index})

    q.join()


run_prague()


# In[ ]:


def get_file_name(crossings, index, type):
    if type == "differential":
        return "knot_" + str(crossings) + "_" + str(index)
    if type == "eig":
        return "knot_" + str(crossings) + "_" + str(index) + "_eigs"
    if type == "laplacian":
        return "knot_" + str(crossings) + "_" + str(index) + "_laspa"


# In[ ]:





# In[ ]:




if __name__ == "__main__":
    import sys
    get_knot_eigs(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))