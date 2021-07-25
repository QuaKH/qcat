#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def parse_pd_code(path, start_index):
    with open(path, "r") as file:
        lines = file.readlines()
        crossing_lines = lines[3::5]
        pd_codes = lines[6::5]
        for index in range(len(pd_codes)):
            pd_code = pd_codes[index][1:].replace("),", ");").strip()
            crossings = crossing_lines[index].split("random_link(")[1].split(",")[0]
            print("compute_pd_code_differential(" + pd_code + ", " + str(crossings) + ", " + str(index+start_index) + ")")

if __name__ == "__main__":
    import sys
    parse_pd_code(sys.argv[1], int(sys.argv[2]))
