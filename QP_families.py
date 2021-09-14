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
    
        QP = QuiverWithPotential(edges, positions=positions)
        QP.add_term_to_potential(edges, coefs, input_format="edges")
        QP.can_mutate[1] = False
        return QP
        
    else:
        m = int(n/2)+1
        edges = [[0, 1],[1, 0],[1, 2],[2, 1],[2, 0],[0, 2]] \
              + [[i, i+1] for i in range(2, m)] \
              + [[i+1, i] for i in range(2, m)] \
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
    
        QP = QuiverWithPotential(edges, potential=[cycles, coefs], positions=positions)
        QP.can_mutate[1] = False
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
