from sage.all import *

from diamond_cutting import *

import time

#L = Lattice([[GF(3)(1), 0], [0, 1]])
#print L
#print L.zero() + L.an_element()

#L_complex = Lattice([[i, 0], [0, 1]])

#L = Lattice([[2, 0], [0, 1]])
#print L([2, 1])

#L = Lattice([[6, 1], [9, 0]])
#print L

#L.voronoi_cell()

cr63 = """
4.47    0.00    0.00    0.00    0.00    0.00    
-3.35    2.96    0.00    0.00    0.00    0.00    
1.34    -3.55    2.37    0.00    0.00    0.00    
-0.22    1.77    -3.55    2.05    0.00    0.00    
-0.22    -0.59    1.77    -3.76    1.53    0.00    
1.34    1.18    0.59    3.07    -2.29    1.32
"""
cr6 = """
 2.4494897427832E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00  
-1.6329931618555E+00   1.8257418583506E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00  
 4.0824829046386E-01  -1.8257418583506E+00   1.5811388300842E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00  
 0.0000000000000E+00   5.4772255750517E-01  -1.8973665961010E+00   1.4491376746189E+00   0.0000000000000E+00   0.0000000000000E+00  
 0.0000000000000E+00   0.0000000000000E+00   6.3245553203368E-01  -1.9321835661586E+00   1.3662601021279E+00   0.0000000000000E+00  
 4.0824829046386E-01   3.6514837167011E-01   3.1622776601684E-01   9.6609178307930E-01  -1.7078251276599E+00   1.3228756555323E+00  
"""
m2a = """
2 0
0 1
"""
m2b = """
-2 -2
2 -2
"""
es6 = """
1.41    0.00    0.00    0.00    0.00    0.00    
-0.71    1.22    0.00    0.00    0.00    0.00    
0.00    -0.82    1.15    0.00    0.00    0.00    
0.00    0.00    -0.87    1.12    0.00    0.00    
0.00    -0.82    -0.58    -0.45    0.37    0.00    
0.00    0.00    -0.87    -0.67    0.55    0.71
"""
lv6 = """
  1.4135770e+000  0.0000000e+000  0.0000000e+000  0.0000000e+000  0.0000000e+000  0.0000000e+000
  3.7281308e-001  1.3635287e+000  0.0000000e+000  0.0000000e+000  0.0000000e+000  0.0000000e+000
 -2.9499631e-001 -2.2516690e-001  1.4037368e+000  0.0000000e+000  0.0000000e+000  0.0000000e+000
 -3.7281308e-001 -2.8456344e-001 -8.7491668e-001  1.0063572e+000  0.0000000e+000  0.0000000e+000
  3.7281308e-001  2.8456344e-001 -1.7307180e-001 -4.5556135e-001  1.2412671e+000  0.0000000e+000
 -7.4569690e-001 -5.6918088e-001  3.4619130e-001 -5.5058489e-001 -6.4855953e-001  6.2011907e-001
 """
matas6 = """
 2.4494897427832E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00  
-4.0824829046386E-01   2.4152294576982E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00  
-4.0824829046386E-01  -4.8304589153965E-01   2.3664319132398E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00  
-4.0824829046386E-01  -4.8304589153965E-01  -5.9160797830996E-01   2.2912878474779E+00   0.0000000000000E+00   0.0000000000000E+00  
-4.0824829046386E-01  -4.8304589153965E-01  -5.9160797830996E-01  -7.6376261582597E-01   2.1602468994693E+00   0.0000000000000E+00  
-4.0824829046386E-01  -4.8304589153965E-01  -5.9160797830996E-01  -7.6376261582597E-01  -1.0801234497346E+00   1.8708286933870E+00  
"""

def test_lattice(data):
    M = []
    for line in data.splitlines():
        row = []
        for value in line.split(' '):
            value = value.strip()
            if value:
                value = QQ(float(value))
                row.append(value)
        if row:
            M.append(row)
    M = matrix(M)
    start = time.clock()
    V = calculate_voronoi_cell(M, radius=3, debug=True)
    stop = time.clock()
    print V.Vrepresentation()
    print V.Hrepresentation()
    print "Computed Voronoi cell in %f seconds" % (stop - start)
    
#test_lattice(m2b)
#test_lattice(cr6)
#test_lattice(lv6)
#test_lattice(matas6)

#test_lattice("1 0 0 0 0\n0 2 1 0 0")

#L = Lattice([[1, 0], [2, 0], [0, 2]])
#print L.voronoi_cell().Vrepresentation()

#L = Lattice(quadratic_form=[[2,0], [0,2]])

"""
#GM = matrix([[0, 3], [3, -1]])
L = Lattice([[1, 0], [0, 1]])
GM = L.basis()
V = calculate_voronoi_cell(matrix(GM))
print V.Hrepresentation()
print V.Vrepresentation()
"""

A = Matrix(RDF, 3, 3, [1.0/n for n in range(1, 10)])
print A.LLL()
