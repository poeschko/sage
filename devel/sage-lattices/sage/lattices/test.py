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
m2 = """
2 0
0 1
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
    V = calculate_voronoi_cell(M, debug=True)
    print V.Vrepresentation()
    
#test_lattice(m2)
test_lattice(cr6) 

"""
#GM = matrix([[0, 3], [3, -1]])
L = Lattice([[1, 0], [0, 1]])
GM = L.basis()
V = calculate_voronoi_cell(matrix(GM))
print V.Hrepresentation()
print V.Vrepresentation()
"""
