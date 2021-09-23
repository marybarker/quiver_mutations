import sys
sys.path.append('../')
from QP_families import *
from mutations import *

# test d2n for a range of n
#a = D2n(6)
#a.draw(time=1)
#a = D2n(8)
#a.draw(time=1)
#a = D2n(10)
#a.draw(time=1)
#a = D2n(12)
#a.draw(time=1)


## cyclic quiver
#a = cyclicQP(4)
#a.draw(time=1)
#a = cyclicQP(5)
#a.draw(time=1)


a = D2n(6)
print(a.mutate(4))
all_mutations = all_mutation_sequences_for_quiver(a)

print(len(all_mutations))
#    print(q.Q1)

print("the quiver is: ")
print(a)
print("and the sequences yielding distinct mutations are: ")
for m in all_mutations:
    print("sequence = " + ", ".join([str(x) for x in m]))
    Q = a.mutate_in_sequence(m, draw=False)
    Q.draw(time=4)
    print("quiver = ", a.mutate_in_sequence(m, draw=False).Q1)

