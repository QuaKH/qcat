#!/usr/bin/env/python
# This script was saved by SnapPy on Sat Sep 11 17:36:26 2021.

crossings=3
# L = random_link(crossings=3, num_components=1, simplify=None)
# Trefoil
L = Link([[1, 4, 2, 5],[3, 6, 4, 1],[5, 2, 6, 3]])
# Unknot
# L = Link([[1,0,0,1]])
L.view(show_crossing_labels=True)
original_PD_code = L.PD_code()
# print initial PD code
print(str(crossings)+";0;"+str(L.PD_code()))
num_iter = 5000
idx = 1
for i in range(num_iter):
    # reset to original random link
    L = Link(original_PD_code)
    # keep track of max number of crossings this knot may have
    max_crossings = crossings
    while max_crossings < 11:
        # perform Reidemeister moves
        L.backtrack(steps=1, prob_type_1=0.33333333333333333333333333333333, prob_type_2=0.33333333333333333333333333333333)
        
        PD_code = L.PD_code()
        print(str(len(PD_code))+";"+str(idx)+";"+str(PD_code))
        
        max_crossings += 2
        idx += 1
