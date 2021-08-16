import sys
sys.path.append('../')
from mutations import *

numRows = 1#3
numCols = 2#4
edges = [[(numCols+1)*r+i, (numCols+1)*r+i+1] for r in range(numRows+1) for i in range(numCols)] \
      + [[(numCols+1)*r+i, (numCols+1)*(r+1)+i] for r in range(numRows) for i in range(numCols+1)] \
      + [[(numCols+1)*r+i, (numCols+1)*(r+1)+i+1] for r in range(numRows) for i in range(numCols)]  # diagonal edges 

hi = 1/numRows
hj = 1/numCols
positions = [[hj*j,hi*i] for i in range(numRows+1) for j in range(numCols+1)]

R = Resolution(edges, vertex_positions=positions)

# visualize the floppable edges of R
for e in range(len(R.edges)):
    color = 'r' if R.can_flop(e) else 'y'
    cs = ['b' for i in range(len(R.edges))]
    cs[e] = color
    kw = {"edge_color":cs}
    R.draw(time=1, **kw)

QP = R.QP()

QP.draw(time=2)
