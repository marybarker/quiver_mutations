"""
This file contains tests for families of QPs that might be useful. 
"""

import sys
sys.path.append('../')
from QP_families import *
from mutations import *

# test d2n for a range of n
a = D2n(6)
a.draw(time=1)
print(a)
a = D2n(8)
a.draw(time=1)
print(a)
a = D2n(10)
a.draw(time=1)
print(a)
a = D2n(12)
a.draw(time=1)
print(a)


# cyclic quiver
a = cyclicQP(4)
a.draw(time=1)
print(a)
a = cyclicQP(5)
a.draw(time=1)
print(a)

a = favorite_example()
a.draw(time=1)
print(a)


a = tetrahedralQP(1)
a.draw(time=1)
print(a)
all_mutations = all_mutation_sequences_for_quiver(a)
print(len(all_mutations))


a = octahedralQP(1)
a.draw(time=1)
print(a)
