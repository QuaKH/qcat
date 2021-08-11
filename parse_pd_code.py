#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def parse_pd_code(path):
    with open(path, "r") as file:
        lines = file.readlines()
        crossing_lines = lines[3::5]
        pd_codes = lines[6::5]
        for index in range(len(pd_codes)):
            pd_code = pd_codes[index][1:].replace("),", ";").replace("(","").replace(")","")
            crossings = crossing_lines[index].split("random_link(")[1].split(",")[0]
            print("compute_pd_code_differential(" + pd_code[:-1] + ", " + str(crossings) + ", " + str(index) + ")")
# parse_pd_code("./pd_code_input.txt")

# def parse_pd_code(path):
#     with open(path, "r") as file:
#         lines = file.readlines()
#         for line in lines:
#             parts = line.split(";")
#             pd_code = parts[2].replace("),", ";").replace("(","").replace(")","")
#             crossings = int(parts[0])
#             index = (crossings - 4) * 50 + int(parts[1])
#             print("compute_pd_code_differential(" + pd_code[:-1] + ", " + str(crossings) + ", " + str(index) + ")")


if __name__ == "__main__":
    import sys
    parse_pd_code(sys.argv[1])

# %%
