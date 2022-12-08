import sys
sys.path.append('../')
from mutations import *

numRows = 3
numCols = 4
edges = [[(numCols+1)*r+i, (numCols+1)*r+i+1] for r in range(numRows+1) for i in range(numCols)] \
      + [[(numCols+1)*r+i, (numCols+1)*(r+1)+i] for r in range(numRows) for i in range(numCols+1)] \
      + [[(numCols+1)*r+i, (numCols+1)*(r+1)+i+1] for r in range(numRows) for i in range(numCols)]  # diagonal edges 

R = Resolution(edges)

## plot a resolution and its flop at an edge
#R.draw()
#R.flop(1).draw()

# pre-set the coordinates of R's vertices so as to make drawing prettier
hi = 1/numRows
hj = 1/numCols
positions = [[hj*j,hi*i] for i in range(numRows+1) for j in range(numCols+1)]
print(positions)

# now animate a sequence of flops (edges to flop specified in list es below)
R.vertex_positions = positions
es = [2,4,7,9,8,4]
#R1 = R.flop_in_sequence(es)

import networkx as nx
import matplotlib.pyplot as plt
#for t in R.triangles:
#    g = nx.Graph()
#    g.add_edges_from(R.edges)
#    vertices = list(zip(*[tuple(positions[i]) for e in t for i in R.edges[e]]))
#    nx.draw_networkx(g, R.vertex_positions)
#    plt.fill(vertices[0], vertices[1])
#    plt.draw()
#    plt.pause(1)
#    plt.clf()
for e in range(len(R.edges)):
    color = 'r' if R.can_flop(e) else 'y'
    cs = ['b' for i in range(len(R.edges))]
    cs[e] = color
    kw = {"edge_color":cs}
    R.draw(time=1, **kw)


# dummy QP to test mutations
#edges = [[0,1],[1,0], # 0,1
#         [0,2],[2,0], # 2,3
#         [0,3],[3,0], # 4,5
#         [1,2],[2,1], # 6,7
#         [1,3],[3,1], # 8,9
#         [2,3],[3,2]] # 10,11
#
#pot = {( 0,6,3):1, ( 0,8,5):-1, 
#       ( 1,4,9):1, ( 1,2,7):-1,
#       (10,5,2):1, (10,9,6):-1,
#       (11,7,8):1, (11,3,4):-1}
#
#QP = QuiverWithPotential(edges, pot)
#print(QP)
#print(QP.mutate(2))
#
#QP = QuiverWithPotential(edges)
#print(QP)
#QP.add_term_to_potential([(0,1,2),(0,1,3),(1,0,3),(1,0,2),(2,3,0),(2,3,1),(3,2,1),(3,2,0)],[1,-1,1,-1,1,-1,1,-1])
#print(QP)
#print(QP.mutate(2))


#edges = [[0,1],[1,2],[2,3],[3,4],[4,5],[5,0],[0,6],[1,6],[2,6],[3,6],[4,6],[5,6]]
#vps = [[0,0],[2,0],[5/3,2/3],[4/3,4/3],[1,2],[.5,1],[1,.5]]
#R = Resolution(edges, vps)
#es = [10,11,11,9,8]
#R1 = R.flop_in_sequence(es)



edges = [[0,1],[1,2],[2,3],[3,4],[4,5],[5,0],[0,6],[1,6],[2,6],[3,6],[4,6],[5,6]]
vps = [[0,0],[2,0],[5/3,2/3],[4/3,4/3],[1,2],[.5,1],[.8,.45]]
R = Resolution(edges, vps)


for e in range(len(R.edges)):
    color = 'r' if R.can_flop(e) else 'y'
    cs = ['b' for i in range(len(R.edges))]
    cs[e] = color
    kw = {"edge_color":cs}
    R.draw(time=1, **kw)


