import sys
sys.path.append('../')
from QP_families import *

a = D2n(6)
a.draw(time=1)

a = D2n(8)
a.draw(time=1)

a = D2n(10)
a.draw(time=1)

a = D2n(12)
a.draw(time=1)


a = cyclicQP(4)
a.draw(time=1)
a = cyclicQP(5)
a.draw(time=1)
