#!/usr/bin/env python3
import csv
import sys
import numpy as np
from scipy.sparse import lil_matrix

################################# Parameters ##################################
if len(sys.argv) < 4:
    print("Usage:",sys.argv[0],"[IMDBFile (in)] [EdgeFile (out)] [DegFile (out)]")
    sys.exit(0)

# IMDB File (input)
IMDBFile = sys.argv[1]
# Edge File (output)
EdgeFile = sys.argv[2]
# Degree File (output)
DegFile = sys.argv[3]

############################# Read the IMDB file ##############################
# [output1]: movie_actor_lst ([movie_id, actor_id]})
# [output2]: movie_id
# [output3]: movie_id
def ReadIMDB():
    # Initialization
    movie_actor_lst = []
    
    # Read the IMDB file
    f = open(IMDBFile, "r")
    for i, line in enumerate(f):
        # Skip the header
        if i < 55:
            continue
        # Read #movies and #actors --> movie_num, actor_num
        elif i == 55:
            lst = line.rstrip("\n").split(" ")
            movie_num = int(lst[0])
            actor_num = int(lst[1])
        # Read the movie-actor list --> movie_actor_lst
        else:
            lst = line.rstrip("\n").split(" ")
            movie_id = int(lst[0])
            actor_id = int(lst[1])
            movie_actor_lst.append([movie_id, actor_id])
    f.close()

    return movie_actor_lst, movie_num, actor_num

#################################### Main #####################################
# Read the IMDB file
movie_actor_lst, movie_num, actor_num = ReadIMDB()

# Make a movie dictionary ({movie_id: [actor_id]}) --> movie_dic
# (Both movie_id and actor_id start with zero)
movie_dic = {}
for i in range(movie_num):
    movie_dic[i] = []
for lst in movie_actor_lst:
    # Both movie_id and actor_id start with zero
    movie_id = lst[0] - 1
    actor_id = lst[1] - 1
    movie_dic[movie_id].append(actor_id)

# Make edges --> edges_lil, edges_csr
print("Making edges.")
edges_lil = lil_matrix((actor_num, actor_num))
deg = np.zeros(actor_num)
for i in range(movie_num):
#for i in range(10):
    if i % 10000 == 0:
        print(i)
    for j in range(len(movie_dic[i])):
        for k in range(j + 1, len(movie_dic[i])):
            # actor indexes --> actor1, actor2
            actor1 = movie_dic[i][j]
            actor2 = movie_dic[i][k]
            if edges_lil[actor1, actor2] == 0:
                edges_lil[actor1, actor2] = 1
#               edges_lil[actor2, actor1] = 1
                deg[actor1] += 1
                deg[actor2] += 1
#edges_csr = edges_lil.tocsr()
a1, a2 = edges_lil.nonzero()
#a1, a2 = edges_csr.nonzero()
print("#edges:", len(a1))

# Output edge information
print("Outputting edge information.")
f = open(EdgeFile, "w")
print("#nodes", file=f)
print(actor_num, file=f)
print("node,node", file=f)
writer = csv.writer(f, lineterminator="\n")
for i in range(len(a1)):
    # actor indexes --> lst
    lst = [a1[i], a2[i]]
    writer.writerow(lst)
f.close()

# Output degree information
print("Outputting degree information.")
f = open(DegFile, "w")
print("node,deg", file=f)
writer = csv.writer(f, lineterminator="\n")
for actor1 in range(actor_num):
    # actor index and her degree --> lst
    lst = [actor1, int(deg[actor1])]
    writer.writerow(lst)
f.close()
