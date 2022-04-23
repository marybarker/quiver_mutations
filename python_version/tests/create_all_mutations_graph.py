import networkx as nx
import sys
sys.path.append('../')
from mutations import *
from QP_families import *

n = 6
a = D2n(n)

all_mutations = [tuple(y) for y in all_mutation_sequences_for_quiver(a)]

parents = []
for i, x in enumerate(all_mutations):
    if len(x) > 0:
        parent[i] = all_mutations.index(x[:-1])
    else:
        parent[i] = -1

g = nx.Graph()
g.add_nodes_from([(i, {"sequence":s}) for i,s in enumerate(all_mutations)])
nodect = len(all_mutations)


for i,p in enumerate(parents):

    if parent[i] > -1:
        g.add_edge([p,i])

    if len(all_mutations[i]) > 1:
        end = list(all_mutations[i])[-2:]
        s1 = list(all_mutations[p])[:-1]+end
        s2 = list(all_mutations[p])[:-1]+end[::-1]

        one = a.mutate_in_sequence(s1, draw=False)
        two = a.mutate_in_sequence(s2, draw=False)

        if str(one.Q1) == str(two.Q1):
            g.add_node((nodect, {"sequence":s2}))
            g.add_edge([i, nodect])
            nodect += 1

