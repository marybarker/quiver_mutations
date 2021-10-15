import sys
import copy
sys.path.append('../')
from mutations import *
from QP_families import *


def try_adding_edges(QP, max_num=3):
    edges = [[x[0], x[1]] for x in QP.Q1]
    p = list(zip(*QP.potential.items()))
    frozen_nodes = QP.frozen_nodes
    
    for i, e in enumerate(QP.Q1):
        for n_e in range(max_num):
            print("adding %d copies of edge %d"%(n_e, i))
            es = [[e[0],e[1]] for j in range(n_e)]

            QP1 = QuiverWithPotential(copy.copy(edges), frozen_nodes=copy.copy(frozen_nodes))
            QP1.add_term_to_potential(p[0], p[1], input_format="edges")
            QP1.add_edges(es)

            all_mutations = list(all_mutation_sequences_for_quiver(QP1))
            print(" -> %d"%len(all_mutations))


def try_perturbing_potential(QP):
    edges = [[x[0], x[1]] for x in QP.Q1]
    p = list(zip(*QP.potential.items()))
    frozen_nodes = QP.frozen_nodes

    all_mutations = list(all_mutation_sequences_for_quiver(QP))
    print("without adding any loops:")
    print(" -> %d"%len(all_mutations))
    for v in QP.Q0:
        if len(QP.loops_at[v]) > 0:

            print("adding a loop at vertex %d to every cycle containing that vertex"%v)
            loop = QP.loops_at[v][0]
            p1  = [list(copy.copy(p[0])), list(copy.copy(p[1]))]

            for it, term in enumerate(p1[0]):
                t = [y[0] for y in QP.arrows_with_head[v] if y[0] in term]

                if len(t) > 0:
                    ti = list(term).index(t[0])
                    p1[0][it] = list(term[:ti]) + [loop[0]] + list(term[ti:])

            QP1 = QuiverWithPotential(copy.copy(edges), frozen_nodes=copy.copy(frozen_nodes))
            QP1.add_term_to_potential([tuple(x) for x in p1[0]], p1[1], input_format="edges")

            all_mutations = list(all_mutation_sequences_for_quiver(QP1))
            print(" -> %d"%len(all_mutations))


n = 6
print("looking at the D_%d case"%(2*n))
a = D2n(n)

try_adding_edges(a)
try_perturbing_potential(a)
