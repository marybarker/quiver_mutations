"""
This file is *supposed* to be for testing out the routine that generates 
all possible triangulations. It's still a work in progress. 
"""
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
