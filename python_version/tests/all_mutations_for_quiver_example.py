"""
This file demos two ways to calculate the distinct quivers that can be obtained by 
mutating at sequences of vertices. 

There are 2 test cases at the moment: 1/6(1,2,3) case and D_12. 
* The 1/6(1,2,3) case is from 'Examples of mutating quivers with potential' by Tom Ducat(pp. 5-6)
* The D_12 case is shown in 'Flops and mutations for Crepant resolutions of polyhedral singularities' By Nolla and Sekiya (p. 12)

"""
import sys
sys.path.append('../')
from mutations import *
from QP_families import *


# Create the 1/6(1,2,3) example
print("looking at the 1/6(1,2,3) case")
QP = favorite_example()
all_mutations = get_all_mutations_from_quiver(QP)
print("There are %d distinct quivers that can be obtained by mutations"%(len(list(all_mutations))))
print("... and they are as follows: ")
for Q in all_mutations:
    Q.draw(time=1)

print("calculating all possible mutations and saving the sequences....")
sequences = all_mutation_sequences_for_quiver(QP)
print("There are %d distinct quivers to be obtained using sequences. This should match previous result"%(len(list(sequences))))

print("... and they are as follows: ")
for s in sequences:
    q = QP.mutate_in_sequence(s,draw=False)
    q.draw(time=1)


n = 6
print("and now looking at the D_%d case"%(2*n))
a = D2n(n)
all_mutations = all_mutation_sequences_for_quiver(a)

print("There are %d distinct quivers for D_%d"%(len(all_mutations), 2*n))
print("and the sequences yielding distinct mutations are: ")
for m in all_mutations:
    print("sequence = " + ", ".join([str(x) for x in m]))
    Q = a.mutate_in_sequence(m, draw=False)
    Q.draw(time=1)
    print("quiver = ", a.mutate_in_sequence(m, draw=False).Q1)

