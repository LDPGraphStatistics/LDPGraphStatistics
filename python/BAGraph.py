#!/usr/bin/env python3
import csv
import sys
import numpy as np
import networkx as nx

################################# Parameters ##################################
if len(sys.argv) < 5:
    print("Usage:",sys.argv[0],"[n (#nodes)] [m (#edges from a new node)] [EdgeFile (out)] [DegFile (out)]")
    sys.exit(0)

# Paramter n (#nodes)
N = int(sys.argv[1])
# Parameter m (#edges from a new node)
M = int(sys.argv[2])
# Edge File (output)
EdgeFile = sys.argv[3]
# Degree File (output)
DegFile = sys.argv[4]

#################################### Main #####################################
# Fix a seed
seed = 1

# Generate a random graph according to the Barabasi-Albert model
G = nx.barabasi_albert_graph(N, M, seed)

# Calculate a degree for each node --> deg
deg = np.zeros(N, dtype=int)
for (i,j) in G.edges():
    deg[i] += 1
    deg[j] += 1

# Output edge information
print("Outputting edge information.")
f = open(EdgeFile, "w")
print("#nodes", file=f)
print(G.number_of_nodes(), file=f)
print("node,node", file=f)
writer = csv.writer(f, lineterminator="\n")
for (i,j) in G.edges():
    lst = [i, j]
    writer.writerow(lst)
f.close()

# Output degree information
print("Outputting degree information.")
f = open(DegFile, "w")
print("node,deg", file=f)
writer = csv.writer(f, lineterminator="\n")
for i in range(N):
    # actor index and her degree --> lst
    lst = [i, deg[i]]
    writer.writerow(lst)
f.close()
