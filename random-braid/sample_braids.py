import random
import snappy

def generate_random_braid_closure(strands, height):
    gens = [x for x in range(1, strands)] + [x for x in range(1 - strands, 0)]
    braid_word = random.choices(gens, k = height)
    return snappy.Link(braid_closure=braid_word)

for strands in range(3,6):
    for height in range(1,12):
        for i in range(300):
            print(generate_random_braid_closure(strands, height).PD_code())