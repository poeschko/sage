from sage.all import *

from diamond_cutting import *

#L = Lattice([[GF(3)(1), 0], [0, 1]])
#print L
#print L.zero() + L.an_element()

#L_complex = Lattice([[i, 0], [0, 1]])

#L = Lattice([[2, 0], [0, 1]])
#print L([2, 1])

#L = Lattice([[6, 1], [9, 0]])
#print L

#L.voronoi_cell()

#GM = matrix([[0, 3], [3, -1]])
L = Lattice([[1, 0], [0, 1]])
GM = L.basis()
V = calculate_voronoi_cell(matrix(GM))
print V.Hrepresentation()
print V.Vrepresentation()