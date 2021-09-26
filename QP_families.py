from mutations import *
import numpy as np

def D2n(n):
    # if n is odd
    if n%2 > 0:
        m = int((n-1)/2)

        edges = [[0, 1],[1, 0],[1, 2],[2, 1],[2, 0],[0, 2]] \
              + [[i, i+1] for i in range(2, m+1)] \
              + [[i+1, i] for i in range(2, m+1)] \
              + [[i, i] for i in range(2, m+1)] \
              + [(m+1,m+1)]
    
        num = m-1
        last = len(edges) - 1
        # have to create potential by  edge index instead of vertices 
        # because of multiloop at vertex m. 
        cycles = [(0, 2, 4), (1, 5, 3)] \
               + [(6+2*num, 2, 3), (6+2*num, 4, 5)] \
               + [(6+2*num+i, 6+num+i-1, 6+i-1) for i in range(1,num)] \
               + [(6+2*num+i, 6+i, 6+num+i) for i in range(num)] \
               + [(last-1, last, last)] \
    
        coefs = [-1,-1] \
              + [1,1] \
              + [1 for i in range(1, num)] \
              + [-1 for i in range(num)] \
              + [-1]
    
        positions = [(0,-1), (0,1)] \
                  + [(x, 0) for x in range(1, m+1)]
    
        QP = QuiverWithPotential(edges, positions=positions, frozen_nodes=[1])
        QP.add_term_to_potential(cycles, coefs, input_format="edges")
        return QP
        
    else:
        m = int(n/2)+1
        edges = [[0, 1],[1, 0],[1, 2],[2, 1],[2, 0],[0, 2]] \
              + [[i, i+1] for i in range(2, m-1)] \
              + [[i+1, i] for i in range(2, m-1)] \
              + [[i, i] for i in range(2, m)] \
              + [[m-1, m],[m, m-1],[m, m+1],[m+1, m],[m+1, m-1],[m-1, m+1]]
    
        cycles = [(0,1,2), (0,2,1), (m-1, m, m+1), (m-1, m+1, m)] \
               + [(2,2,1), (2,2,0), (m-1,m-1,m), (m-1,m-1,m+1)] \
               + [(i,i,i-1) for i in range(3, m)] \
               + [(i,i,i+1) for i in range(2, m-1)]
    
        coefs = [-1,-1,-1,-1] \
              + [1,1,-1,-1] \
              + [1 for i in range(3, m)] \
              + [-1 for i in range(2, m-1)]
    
        positions = [(0,-1), (0,1)] \
                  + [(x, 0) for x in range(1, m-1)] \
                  + [(m-1, 1), (m-1, -1)]
    
        QP = QuiverWithPotential(edges, potential=[cycles, coefs], positions=positions, frozen_nodes=[1])
        return QP


def cyclicQP(n):
    edges = [[i, (i+1)%n] for i in range(n)] \
          + [[i, (i+n-1)%n] for i in range(n)] \
          + [[i, i] for i in range(n)]


    cycles = [(i,(i+1)%n, i) for i in range(n)] \
           + [(i,(i+n-1)%n, i) for i in range(n)]
    coefs = [1 for i in range(n)] + [-1 for i in range(n)]

    angle = 2*np.pi / n

    positions = [(np.sin(angle*i), np.cos(angle*i)) for i in range(n)]

    return QuiverWithPotential(edges, potential=[cycles, coefs], positions=positions)


def favorite_example():
    # this is the 1/6(1,2,3) singularity quiver. It's used so frequently I don't want to keep typing it. 
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
    positions = [[-1,0], [0,1], [1,0], [1,-1], [0,-2], [-1,-1]]

    return QuiverWithPotential(edges, [cycles, coefs], positions=positions)


def tetrahedralQP(w=1):
    edges = [[0,3],[3,0],[1,3],[3,1],[2,3],[3,2],[3,3],[3,3]]
    cycles = [(6,1,0), (6,3,2), (6,5,4), (6,6,6), (7,1,0), (7,3,2), (7,5,4), (7,7,7)]
    coefs = [1,w,w**2,-1/3,-1,-w**2,-w,1/3]
    positions = [[0,-2],[-1,1],[1,1],[0,0]]

    QP = QuiverWithPotential(edges, positions=positions, frozen_nodes=[0])
    QP.add_term_to_potential(cycles, coefs, input_format="edge_order")
    return QP


def octahedralQP(w):
    edges = [[0,3],[3,0],
             [3,2],[2,3],
             [2,4],[4,2],
             [3,4],[4,3],
             [1,4],[4,1],
             [3,3],[4,4]]
    cycles = [(3,3,0),(3,3,2),(3,3,4),(3,3,3),(4,4,1),(4,4,2),(4,4,3),(4,4,4),(3,4,2),(4,3,2)]
    coefs = [1,-1,-1,-1/3,1,-1,-1,1/3,w**2-w,w**2-w]
    positions = [[-2,0],[2,0],[0,2],[-1,0],[1,0]]
    return QuiverWithPotential(edges, [cycles,coefs], positions=positions)
