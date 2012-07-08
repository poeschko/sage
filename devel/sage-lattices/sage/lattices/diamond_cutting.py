from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import matrix, identity_matrix
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
    
    # coerce to floats 
    GM = GM.N()    
    C = float(C)
    if debug:
        print "Cut\n%s\nwith radius %s" % (GM, C)
    
    dim = GM.dimensions()
    assert dim[0] == dim[1]
    dim = dim[0]
    T = [0] * dim
    U = [0] * dim
    x = [0] * dim
    L = [0] * dim
    
    # calculate the Gram matrix
    q = matrix([[sum(GM[i][k] * GM[j][k] for k in range(dim)) for j in range(dim)] for i in range(dim)])
    if debug:
        print "q:\n%s" % q.N()
    # apply Cholesky/Jacobi decomposition
    q = jacobi(q)
    if debug:
        print "q:\n%s" % q.N()
    
    i = dim - 1
    T[i] = C
    U[i] = 0
    
    new_dimension = True
    cut_count = 0
    inequalities = []
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
        elif i > 0:
            T[i - 1] = T[i] - q[i][i] * (x[i] + U[i]) ** 2
            i -= 1
            U[i] = 0
            for j in range(i + 1, dim):
                U[i] += q[i][j] * x[j]
            new_dimension = True
        else:
            if all(elmt == 0 for elmt in x):
                break
            hv = [0] * dim
            for k in range(dim):
                for j in range(dim):
                    hv[k] += x[j] * GM[j][k]
            hv = vector(hv)
                    
            for hv in [hv, -hv]:
                cut_count += 1
                if debug:
                    print "\n%d) Cut using normal vector %s" % (cut_count, hv)
                hv = [QQ(round(elmt, 6)) for elmt in hv]
                inequalities.append(plane_inequality(hv))
                #cut = Polyhedron(ieqs=[plane_inequality(hv)])
                #V = V.intersection(cut)
           
    if debug:
        print "Final cut"     
    cut = Polyhedron(ieqs=inequalities)
    V = V.intersection(cut)
    
    if debug:
        print "End"
    
    return V
    
def calculate_voronoi_cell(basis, radius=None, debug=False):
    " Calculate the Voronoi cell of the lattice defined by basis "
    
    dim = basis.dimensions()
    artificial_length = None
    if dim[0] < dim[1]:
        # introduce "artificial" basis points (representing infinity)
        artificial_length = ceil(max(abs(v) for v in basis)) * 2
        additional_vectors = identity_matrix(dim[1]) * artificial_length
        basis = basis.stack(additional_vectors)
        # LLL-reduce to get quadratic matrix
        basis = basis.LLL()
        basis = matrix([v for v in basis if v])
        dim = basis.dimensions()
    assert dim[0] == dim[1]
    basis = basis / 2
    
    ieqs = []
    for v in basis:
        ieqs.append(plane_inequality(v))
        ieqs.append(plane_inequality(-v))
    Q = Polyhedron(ieqs=ieqs)
    
    # twice the length of longest vertex in Q is a safe choice
    if radius is None:
        radius = 2 * max(abs(v.vector()).N() for v in Q.vertex_generator())
    
    V = diamond_cut(Q, basis, radius, debug=debug)
    
    if artificial_length is not None:
        # remove inequalities introduced by artificial basis points
        H = V.Hrepresentation()
        H = [v for v in H if all(not V._is_zero(v.A() * w / 2 - v.b() and
            not V._is_zero(v.A() * (-w) / 2 - v.b())) for w in additional_vectors)]
        V = Polyhedron(ieqs=H)
        
    return V