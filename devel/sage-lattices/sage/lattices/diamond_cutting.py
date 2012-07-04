from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.rational_field import QQ
#from sage.rings.infinity import PlusInfinity

from math import sqrt, floor, ceil

def plane_inequality(v):
    """ Return the inequality for points on the same side as the origin
    with respect to the plane through v normal to v.
    
    sage: from sage.lattices.diamond_cutting import plane_inequality
    sage: ieq = plane_inequality([1, -1]); ieq
    [2, -1, 1]
    sage: ieq[0] + vector(ieq[1:]) * vector([1, -1])
    0
    """
    
    v = vector(v)
    c = -v * v
    if c < 0:
        c, v = -c, -v
    return [c] + list(v)

def jacobi(M):
    " Cholesky/Jacobi decomposition of M "
    
    dim = M.dimensions()
    assert dim[0] == dim[1]
    dim = dim[0]
    q = [list(row) for row in M]
    for i in range(dim - 1):
        for j in range(i + 1, dim):
            q[j][i] = q[i][j]
            q[i][j] = q[i][j] / q[i][i]
        for k in range(i + 1, dim):
            for l in range(k, dim):
                q[k][l] = q[k][l] - q[k][i] * q[i][l]
    for i in range(1, dim):
        for j in range(i):
            q[i][j] = 0
    return matrix(q)

def diamond_cut(V, GM, C, debug=False):
    " Perform diamond cutting on Polyhedron V with basis matrix GM and radius C "
    
    C = float(C)
    if debug:
        print "Cut\n%s\nwith radius %s" % (GM, C)
    
    dim = len(GM[0])
    T = [0] * dim
    U = [0] * dim
    x = [0] * dim
    L = [0] * dim
    
    #q = GM * GM
    q = matrix([[sum(GM[i][k] * GM[j][k] for k in range(dim)) for j in range(dim)] for i in range(dim)])
    if debug:
        print "q:\n%s" % q.N()
    #q = q.cholesky_decomposition()
    q = jacobi(q)
    if debug:
        print "q:\n%s" % q.N()
    
    i = dim - 1
    T[i] = C
    U[i] = 0
    
    new_dimension = True
    cut_count = 0
    while True:
        if debug:
            print "Dimension: %d" % i
        if new_dimension:
            Z = sqrt(T[i] / q[i][i])
            if debug:
                print "Z: %s" % Z
            L[i] = int(floor(Z - U[i]))
            if debug:
                print "L: %s" % L
            x[i] = int(ceil(-Z - U[i]) - 1)
            new_dimension = False
    
        x[i] += 1
        if debug:
            print "x: %s" % x
        if x[i] > L[i]:
            i += 1
            continue
        if i > 0:
            T[i - 1] = T[i] - q[i][i] * (x[i] + U[i]) ** 2
            i -= 1
            U[i] = 0
            for j in range(i + 1, dim):
                U[i] += q[i][j] * x[j]
            new_dimension = True
            continue
        #else:
        if all(elmt == 0 for elmt in x):
            break
        hv = [0] * dim
        for k in range(dim):
            for j in range(dim):
                hv[k] += x[j] * GM[j][k].N()
        hv = vector(hv)
                
        for hv in [hv, -hv]:
            cut_count += 1
            if debug:
                print "\n%d) Cut using normal vector %s" % (cut_count, hv)
            hv = [QQ(elmt) for elmt in hv]
            cut = Polyhedron(ieqs=[plane_inequality(hv)])
            V = V.intersection(cut)
        
    if debug:
        print "End"
    
    return V
    
def calculate_voronoi_cell(basis, debug=False):
    " Calculate the Voronoi cell of the lattice defined by basis "
    
    basis = basis / 2
    
    ieqs = []
    for v in basis:
        ieqs.append(plane_inequality(v))
        ieqs.append(plane_inequality(-v))
    Q = Polyhedron(ieqs=ieqs)
    
    # twice the length of longest vertex in Q is a safe choice
    """radius = 0
    for v in Q.vertex_generator():
        #print v
        l = abs(v.vector()).N()
        if l > radius:
            radius = l
        #print "%s: %d" % (v, l)
        #if l < 100:
    radius *= 2"""
            
    #radius = 6
    radius = 2 * max(abs(v.vector()).N() for v in Q.vertex_generator())
    #radius = max(abs(vector(v)) for v in basis)#
    #radius = 6
    
    V = diamond_cut(Q, basis, radius, debug=debug)
    return V
