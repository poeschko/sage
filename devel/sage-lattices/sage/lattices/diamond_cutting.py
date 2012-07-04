from sage.geometry.polyhedron.constructor import Polyhedron
from sage.modules.free_module_element import vector
from sage.rings.rational_field import QQ

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

def diamond_cut(V, GM, C, debug=False):
    if debug:
        print "Cut\n%s\nwith radius %s" % (GM, C)
    
    dim = len(GM[0])
    T = [0] * dim
    U = [0] * dim
    x = [0] * dim
    L = [0] * dim
    
    #print GM
    #print GM * GM
    
    q = (GM * GM) #.cholesky_decomposition()
    
    i = dim - 1
    T[i] = C
    U[i] = 0
    
    #def process_dimension(i):
    new_dimension = True
    cut_count = 0
    while True:
        if debug:
            print "Dimension: %d" % i
        if new_dimension:
            Z = sqrt(T[i] / q[i][i])
            if debug:
                print "Z: %f" % Z
            L[i] = floor(Z - U[i])
            if debug:
                print "L: %s" % L
            x[i] = ceil(-Z - U[i]) - 1
            new_dimension = False
    
        x[i] += 1
        if debug:
            print "x: %s" % x
        if x[i] > L[i]:
            i += 1
        elif i > 0:
            sum = x[i] + U[i]
            T[i - 1] = T[i] - q[i][i] * (sum ** 2)
            i -= 1
            U[i] = 0
            for j in range(i, dim):
                U[i] += q[i][j] * x[j]
            new_dimension = True
        else:
            #terminate = True
            #for k in range(dim):
            #    terminate = terminate and x[k] == 0
            #if terminate:
            if all(elmt == 0 for elmt in x):
                break
            hv = [0] * dim
            for k in range(dim):
                for j in range(dim):
                    hv[k] += x[j] * GM[j][k]
            cut_count += 1
            if debug:
                print "\n%d) Cut using normal vector %s" % (cut_count, hv)
            hv = [QQ(elmt) for elmt in hv]
            cut = Polyhedron(ieqs=[plane_inequality(hv)])
            V = V.intersection(cut)
    if debug:
        print "End"
    
    return V
    
def calculate_voronoi_cell(basis):
    ieqs = []
    for v in basis:
        ieqs.append(plane_inequality(v / 2))
        ieqs.append(plane_inequality(-v / 2))
    Q = Polyhedron(ieqs=ieqs)
    #print Q
    
    # twice the length of longest vertex in Q is a safe choice
    radius = 2 * max(abs(v.vector()) for v in Q.vertex_generator())
    
    V = diamond_cut(Q, basis, radius)
    return V
