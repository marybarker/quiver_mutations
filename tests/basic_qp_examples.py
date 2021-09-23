import sys
sys.path.append('../')
from mutations import *

# dummy QP to test mutations
edges = [[0,1],[1,0], # 0,1
         [0,2],[2,0], # 2,3
         [0,3],[3,0], # 4,5
         [1,2],[2,1], # 6,7
         [1,3],[3,1], # 8,9
         [2,3],[3,2]] # 10,11

# create cycles, where we interpret (0,1,2) as the cycle x01->x12->x20, etc.
cycles = [(0,1,2),(0,1,3),(1,0,3),(1,0,2),(2,3,0),(2,3,1),(3,2,1),(3,2,0)]

# these are the coefficients for the cycles defined above
coefs = [1,-1,1,-1,1,-1,1,-1]

QP = QuiverWithPotential(edges, [cycles, coefs])

print(QP)
print(QP.mutate(2))

# draw without specifying vertex positions
QP.draw()

# now set vertex positions and re-draw
QP.positions = [[-1,-1], [-1,1], [1,1], [1,-1]]
QP.mutate(2).draw()

QP.mutate_in_sequence([2])
