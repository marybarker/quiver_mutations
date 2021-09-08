import sys
sys.path.append('../')
from mutations import *

edges = [[0,1],[1,0],
         [1,2],[2,1],
         [0,2],[2,0], 
         [2,3],[3,2], 
         [3,3],
         [3,4],[4,3],
         [4,4],
         [4,5],[5,4],
         [5,5],[5,5],
         [5,0],[0,5]]

cycles = [(0,1,2),(0,2,1),(3,3,4),(1,2,3,2),(0,5,0,1),(2,3,3)]
coefs = [1,1,1,1,1,1]

QP = QuiverWithPotential(edges, [cycles, coefs])

# set vertex positions
QP.positions = [[-1,0], [0,1], [1,0], [1,-1], [0,-2], [-1,-1]]

all_mutations = get_all_mutations_from_quiver(QP)
for Q in all_mutations:
    Q.draw(time=1)

print("and now we're doing it the other way!")
sequences = all_mutation_sequences_for_quiver(QP)
for s in sequences:
    q = QP.mutate_in_sequence(s,draw=False)
    q.draw(time=1)
