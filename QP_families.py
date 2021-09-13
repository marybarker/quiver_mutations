from mutations import *

def D2n(n):

    # if n is odd
    if n%4 > 0:
        m = int((n-1)/2)
        print("m is odd", m)
        vertices = list(range(m+2))

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
        return QP
        
    else:
        m = int(n/2)
        print("m is even", m)
        vertices = list(range(m+2))
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
    

        return QuiverWithPotential(edges, potential=[cycles, coefs], positions=positions)
