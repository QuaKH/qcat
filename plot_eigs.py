import matplotlib.pyplot as plt
import math

def plot_knot_eigs(crossings, index):
    path = "./eigs/knot_" + str(crossings) + "_" + str(index) + "_eigs"
    dict = {}
    i_vals = []
    j_vals = []
    ratios = []

    with open(path) as reader:
        line = reader.readline().strip().split(" ")
        while line != ['']:
            i_vals.append(int(line[0]))
            j_vals.append(int(line[1]))

            eigs = []
            for eig_val in line[2:]:
                eigs.append(float(eig_val))
            
            min_idx = 0
            while abs(eigs[min_idx]) < 1e-4: # math.isclose(0,, rel_tol=1e-4):
                min_idx += 1

            ratios.append(eigs[min_idx] / eigs[-1])

            line = reader.readline().strip().split(" ")
    return i_vals, j_vals, ratios

def get_min_eig_ratio(crossings):
    ratios = []

    for index in range(0 + (crossings - 4) * 50, 50 + (crossings - 4) * 50):
        i_vals, j_vals, ratios = plot_knot_eigs(crossings, index)
    return min(ratios)

if __name__ == "__main__":
    max_ratios = []
    for crossings in range(4,13):
        max_ratios.append(get_min_eig_ratio(crossings))

    plt.scatter(range(4,13), max_ratios)
    plt.savefig("min_link_ratios.png")