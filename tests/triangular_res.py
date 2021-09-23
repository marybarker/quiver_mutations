"""
This file creates the 1/6(1,2,3) resolution (modulo I've got vertex positions off)
and runs through which edges can be flopped, before performing the sequence of 
flops as in page 3 of Tom Ducat's 'Examples of mutating quivers with potential'
"""
import sys
sys.path.append('../')
from mutations import *

edges = [[0,1],[1,2],[2,3],[3,4],[4,5],[5,0],[0,6],[1,6],[2,6],[3,6],[4,6],[5,6]]
vps = [[0,0],[2,0],[5/3,2/3],[4/3,4/3],[1,2],[.5,1],[.8,.45]]
R = Resolution(edges, vps)


for e in range(len(R.edges)):
    color = 'r' if R.can_flop(e) else 'y'
    cs = ['b' for i in range(len(R.edges))]
    cs[e] = color
    kw = {"edge_color":cs}
    R.draw(time=1, **kw)


es = [10,11,11,9,8]
R1 = R.flop_in_sequence(es)

QP = R.QP()

QP.draw(time=2)
