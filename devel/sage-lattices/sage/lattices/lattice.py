"""
Lattices

This module provides the base class :class:`Lattice_with_basis` from which
all lattices in Sage derive, as well as a selection of more
specific base classes.

The class inheritance hierarchy is:

- :class:`Lattice_with_basis`

  - :class:`Lattice_ZZ`

    - :class:`Lattice_ZZ_in_RR`
    
Real-valued vectors are coerced to rationals when used in lattices
(see discussion at
https://groups.google.com/forum/?fromgroups#!topic/sage-devel/omE8nI2HFbg).
Complex-valued vectors do not work in lattices.

AUTHORS:

- Jan Poeschko (2012-05-26): initial version

"""

import sage.matrix.matrix_space

from sage.misc.latex import latex

from sage.modules.free_module import FreeModule_submodule_with_basis_pid, element_class
from sage.modules.free_module_element import vector
from sage.modules.free_quadratic_module import FreeQuadraticModule_ambient_pid
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.libs.pari.gen import pari
from sage.matrix.matrix_space import MatrixSpace
from sage.symbolic.ring import SymbolicRing
        
from diamond_cutting import calculate_voronoi_cell

class Lattice_with_basis(FreeQuadraticModule_ambient_pid):
#(FreeModule_submodule_with_basis_pid):
    """
    Construct a general lattice over a PID with a given basis.
    
    INPUT:
    
    - ``ambient`` -- ambient free module over a principal ideal domain `R`,
      i.e. `R^n`;

    - ``basis`` -- list of elements of `K^n`, where `K` is the fraction field
      of `R`.
      
    EXAMPLES::
    
        sage: Lattice([[2, 0], [0, 1]], base_ring=Zp(5))
        Lattice of degree 2 and rank 2 over 5-adic Ring with capped relative precision 20
        Inner product matrix:
        [4 0]
        [0 1]
        Basis matrix:
        [2 0]
        [0 1]
        
    Real-valued (induced) inner product matrix::
    
        sage: L = Lattice([[1.0, 0, 0], [0, 1, 0]]); L
        Lattice of degree 3 and rank 2 over Integer Ring
        Inner product matrix:
        [ 1.00000000000000 0.000000000000000]
        [0.000000000000000  1.00000000000000]
        Basis matrix:
        [ 1.00000000000000 0.000000000000000 0.000000000000000]
        [0.000000000000000  1.00000000000000 0.000000000000000]
    """
    
    def __init__(self, base_ring, rank, inner_product_matrix, embedded_basis):
        """
        See :class:`Lattice_with_basis` for documentation.
        
        TESTS::
        
            sage: L = Lattice([[1, 0], [0, 1]], base_ring=Zp(5))
            sage: L.random_element() # random
            (-1, 2)
            sage: TestSuite(L).run()
        """
        super(Lattice_with_basis, self).__init__(base_ring, rank, inner_product_matrix)
        #super(Lattice_with_basis, self).__init__(ambient, basis, inner_product_matrix)
        #super(Lattice_with_basis, self).__init__(ambient, basis, echelonize=False,
        #    already_echelonized=True)
        self._embedded_basis = matrix(embedded_basis).rows()
        if not self._embedded_basis:
            raise ValueError("basis must not be empty")
        
    def degree(self):
        return len(self._embedded_basis[0])
        
    def embedded_basis(self):
        return self._embedded_basis
        
    def embedded_basis_matrix(self):
        return matrix(self._embedded_basis)
    
    def to_embedded(self, vector):
        return sum(v * b for v, b in zip(vector, self.embedded_basis()))
    
    def determinant(self):
        return self.inner_product_matrix().determinant()
    
    def discriminant(self):
        return abs(self.determinant())
        
    def _repr_(self):
        """
        Text representation of this lattice.
        
        TESTS::
        
            sage: Lattice([[2, 0], [0, 1]], base_ring=Zp(5))
            Lattice of degree 2 and rank 2 over 5-adic Ring with capped relative precision 20
            Inner product matrix:
            [4 0]
            [0 1]
            Basis matrix:
            [2 0]
            [0 1]        
            sage: Lattice([[1.0, 0, 0], [0, 1, 0]])
            Lattice of degree 3 and rank 2 over Integer Ring
            Inner product matrix:
            [ 1.00000000000000 0.000000000000000]
            [0.000000000000000  1.00000000000000]
            Basis matrix:
            [ 1.00000000000000 0.000000000000000 0.000000000000000]
            [0.000000000000000  1.00000000000000 0.000000000000000]
        """
        return "Lattice of degree %s and rank %s over %s\nInner product matrix:\n%s\nBasis matrix:\n%s"%(
            self.degree(), self.rank(), self.base_ring(), self.inner_product_matrix(), self.embedded_basis_matrix())
        
class Lattice_ZZ(Lattice_with_basis):
    """
    Construct a ZZ-lattice that is embedded in a real-valued space
    (e.g. RR^n, QQ^n, or ZZ^n).
    
    EXAMPLES::
    
        sage: L = Lattice([[1, 0, 0], [0, 1, 0]]); L
        ZZ-lattice of degree 3 and rank 2
        Inner product matrix:
        [1 0]
        [0 1]
        Basis matrix:
        [1 0 0]
        [0 1 0]
        
    By default, the basis is reduced using the LLL algorithm::
    
        sage: L = Lattice([[6, 1], [9, 0]]); L
        ZZ-lattice of degree 2 and rank 2
        Inner product matrix:
        [ 9 -3]
        [-3 10]
        Basis matrix:
        [ 0  3]
        [ 3 -1]
        sage: L.discriminant()
        81
        
    However, you can prevent this::
    
        sage: Lattice([[6, 1], [9, 0]], reduce=False)
        ZZ-lattice of degree 2 and rank 2
        Inner product matrix:
        [37 54]
        [54 81]
        Basis matrix:
        [6 1]
        [9 0]
        sage: L.discriminant()
        81
    """
    def __init__(self, base_ring, rank, inner_product_matrix, embedded_basis, reduce=True):
        """
        See :class:`Lattice_ZZ` for documentation.
        
        TESTS::
        
            sage: L = Lattice([[2, 1, 0], [0, 1, 0]])
            sage: TestSuite(L).run()
        """
        '''if reduce_basis:
            # coerce vectors to QQ^n (QQ is the fraction field of ZZ)
            # so that LLL algorithm can be applied
            degree = len(basis[0])
            fraction_field = QQ ** degree
            basis = [fraction_field(x) for x in basis]
            
            # create basis matrix and apply LLL algorithm 
            basis_matrix = matrix(basis)
            basis_matrix = basis_matrix.LLL()
            
            # remove 0-vectors from basis
            basis = list(v for v in basis_matrix if v)
            self.__reduced_basis = basis
        else:
            self.__reduced_basis = None'''
        #super(Lattice_ZZ, self).__init__(basis)
        self.__reduced_inner_product_matrix = None
        self.__reduced_embedded_basis_matrix = None
        if reduce:
            inner_product_matrix, embedded_basis_matrix = self.__reduce(inner_product_matrix,
                matrix(embedded_basis))
            embedded_basis = embedded_basis_matrix.rows()
        
        super(Lattice_ZZ, self).__init__(base_ring, rank, inner_product_matrix, embedded_basis)
        
        self.__voronoi_cell = None  # cached result of voronoi_cell
        #self.__reduced_basis = None
        self.__shortest_vectors = None
        
    '''def embedded_basis(field = None):
        L = self.gram_matrix()
        if field:
            L = L.change_ring(field)
        return L.cholesky().rows()'''
        
    def __reduce(self, inner_product_matrix, embedded_basis_matrix):
        U = inner_product_matrix.LLL_gram()
        reduced_ipm = U.transpose() * inner_product_matrix * U
        reduced_ebm = U.transpose() * embedded_basis_matrix
        self.__reduced_inner_product_matrix = reduced_ipm
        self.__reduced_embedded_basis_matrix = reduced_ebm
        return reduced_ipm, reduced_ebm
        
    def reduced_embedded_basis_matrix(self):
        """
        Return an LLL-reduced basis for this lattice.
        
        EXAMPLES::
        
            sage: L = Lattice([[6, 1], [9, 0]], reduce=False); L
            ZZ-lattice of degree 2 and rank 2
            Inner product matrix:
            [37 54]
            [54 81]
            Basis matrix:
            [6 1]
            [9 0]
            sage: L.reduced_embedded_basis()
            [(0, 3), (3, -1)]
        """
        if self.__reduced_embedded_basis_matrix is not None:
            return self.__reduced_embedded_basis_matrix
        ipm, ebm = self.__reduce(self.inner_product_matrix(), self.embedded_basis_matrix())
        return ebm
        """if self.__reduced_basis is None:
            basis_matrix = self.basis_matrix()
            basis_matrix = basis_matrix.LLL()
            basis = list(v for v in basis_matrix if v)
            self.__reduced_basis = basis
        return self.__reduced_basis"""
        
    def reduced_embedded_basis(self):
        return self.reduced_embedded_basis_matrix().rows()
        
    def _repr_(self):
        """
        Text representation of this lattice.
        
        TESTS::
        
            sage: Lattice([[1, 0, 0], [0, 1, 0]])
            ZZ-lattice of degree 3 and rank 2
            Inner product matrix:
            [1 0]
            [0 1]
            Basis matrix:
            [1 0 0]
            [0 1 0]
        """
        return "ZZ-lattice of degree %s and rank %s\nInner product matrix:\n%s\nBasis matrix:\n%s"%(
            self.degree(), self.rank(), self.inner_product_matrix(), self.embedded_basis_matrix())
        
    def _shortest_vectors(self, max_count=None, max_length=0):
        default = max_length == 0 and max_count == 0
        if default and self.__shortest_vectors is not None:
            return self.__shortest_vectors
        qf = self.gram_matrix()
        #if max_length is None:
        #    # choose trivial upport bound for vector length
        #    max_length = sum(sum(x) ** 2 for x in self.embedded_basis())
        if max_count is None:
            max_count = self._shortest_vectors(max_length=max_length, max_count=0)[0] #self.degree()
        #if self.base_ring() == QQ:
        #    flag = 2
        #else:
        #    flag = 0
        #flag = 2 # allow non-integral entries 
        count, length, vectors = pari(qf).qfminim(0, max_count)
        #vectors = vectors.python()
        #return [self.linear_combination_of_basis(v) for v in vectors.columns()]
        result = count, length, vectors.python().columns()
        if default:
            self.__shortest_vectors = result
        return result
    
    def shortest_vectors_count(self, max_length=0):
        """
        Find the number of shortest vectors in the lattice.
        """
        count, length, vectors = self._shortest_vectors(max_count=0, max_length=max_length)
        return count
    
    def shortest_vectors_length(self):
        """
        Find the length of the shortest vectors in the lattice.
        """
        count, length, vectors = self._shortest_vectors(max_count=0)
        return length
    
    def shortest_vectors(self, max_count=None, max_length=0):
        """
        Find shortest vectors using Pari's Fincke-Pohst algorithm.
        
        INPUT:
          
        - ``max_count`` - limit on the number of vectors returned
          (default: no limit);
        
        - ``max_length`` - maximum length to consider in search
          (default of 0 means search for shortest vectors).
          
        OUTPUT:
        
        A list of at most ``max_count`` shortest vectors.
        Vectors are given in their integer representations with respect to the
        embedded basis.
        
        EXAMPLES::
        
            sage: L = Lattice([[2, 0], [0, 3]])
            sage: L.shortest_vectors()
            [(1, 0)]
            sage: map(L.to_embedded, L.shortest_vectors())
            [(2, 0)]
            
        Note that the given basis might be reduced, leading to unexpected results::
        
            sage: L = Lattice([[2, 1], [1, 1]])
            sage: L.shortest_vectors()
            [(0, 1), (1, 0)]
            sage: map(L.to_embedded, L.shortest_vectors())
            [(0, -1), (-1, 0)]
        """
        count, length, vectors = self._shortest_vectors(max_count=max_count,
            max_length=max_length)
        return vectors
    
    def voronoi_cell(self, radius=None):
        """
        Compute the Voronoi cell of a lattice, returning a Polyhedron.
        
        INPUT:
        
        - ``radius`` -- radius of ball containing considered vertices
          (default: automatic determination).
          
        OUTPUT:
        
        The Voronoi cell as a Polyhedron instance.
        
        The result is cached so that subsequent calls to this function
        return instantly.
        
        EXAMPLES::
        
            sage: L = Lattice([[1, 0], [0, 1]])
            sage: V = L.voronoi_cell()
            sage: V.Vrepresentation()
            (A vertex at (1/2, -1/2), A vertex at (1/2, 1/2), A vertex at (-1/2, 1/2), A vertex at (-1/2, -1/2))
            
        Lattices not having full dimension are handled as well::
        
            sage: L = Lattice([[2, 0, 0], [0, 2, 0]])
            sage: V = L.voronoi_cell()
            sage: V.Hrepresentation()
            (An inequality (-1, 0, 0) x + 1 >= 0, An inequality (0, -1, 0) x + 1 >= 0, An inequality (1, 0, 0) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0)
        
        "Over-dimensional" lattices are reduced first::
        
            sage: L = Lattice([[1, 0], [2, 0], [0, 2]])
            sage: L.voronoi_cell().Vrepresentation()
            (A vertex at (1/2, -1), A vertex at (1/2, 1), A vertex at (-1/2, 1), A vertex at (-1/2, -1))
            
        ALGORITHM:
        
        Uses parts of the algorithm from [Vit1996].
        
        REFERENCES:
        
        .. [Vit1996] E. Viterbo, E. Biglieri. Computing the Voronoi Cell
          of a Lattice: The Diamond-Cutting Algorithm.
          IEEE Transactions on Information Theory, 1996.
        """
        if self.__voronoi_cell is None:
            basis_matrix = self.reduced_embedded_basis_matrix()
            self.__voronoi_cell = calculate_voronoi_cell(basis_matrix, radius=radius)
        return self.__voronoi_cell
    
    def voronoi_relevant_vectors(self):
        """
        Compute the vectors inducing the Voronoi cell.
        
        OUTPUT:
        
        The list of Voronoi relevant vectors.
        
        EXAMPLES::
        
            sage: L = Lattice([[3, 0], [4, 0]])
            sage: L.voronoi_relevant_vectors()
            [(-1, 0), (1, 0)]
        """
        V = self.voronoi_cell()
        
        def defining_point(ieq):
            """
            Compute the point defining an inequality.
            
            INPUT:
            
            - ``ieq`` - an inequality in the form [c, a1, a2, ...]
              meaning a1 * x1 + a2 * x2 + ... <= c
              
            OUTPUT:
            
            The point orthogonal to the hyperplane defined by ``ieq``
            in twice the distance from the origin.
            """
            c = ieq[0]
            a = ieq[1:]
            n = sum(y ** 2 for y in a)
            return vector([2 * y * c / n for y in a])

        return [defining_point(ieq) for ieq in V.inequality_generator()]
    
    def closest_vector(self, t):
        """
        Compute the closest vector in the lattice to a given vector.
        
        INPUT:
        
        - ``t`` -- the target vector to compute the closest vector to.
        
        OUTPUT:
        
        The vector in the lattice closest to ``t``.
        
        EXAMPLES::
        
            sage: L = Lattice([[1, 0], [0, 1]])
            sage: L.closest_vector((-6, 5/3))
            (-6, 2)
            
        ALGORITHM:
        
        Uses the algorithm from [Mic2010].
        
        REFERENCES:
        
        .. [Mic2010] D. Micciancio, P. Voulgaris. A Deterministic Single
          Exponential Time Algorithm for Most Lattice Problems based on
          Voronoi Cell Computations.
          Proceedings of the 42nd ACM Symposium Theory of Computation, 2010.
        """
        voronoi_cell = self.voronoi_cell()
        
        def projection(M, v):
            Mt = M.transpose()
            P = Mt * (M * Mt) ** (-1) * M
            return P * v
        
        t = projection(matrix(self.basis()), vector(t))
        
        def CVPP_2V(t, V, voronoi_cell):
            t_new = t
            while not voronoi_cell.contains(t_new.list()):
                v = max(V, key=lambda v: t_new * v / v.norm() ** 2)
                t_new = t_new - v
            return t - t_new
            
        V = self.voronoi_relevant_vectors()
        t = vector(t)
        p = 0
        while not (ZZ(2 ** p) * voronoi_cell).contains(t):
            p += 1
        t_new = t
        i = p
        while i >= 1:
            V_scaled = [v * (2 ** (i - 1)) for v in V]
            t_new = t_new - CVPP_2V(t_new, V_scaled, ZZ(2 ** (i - 1)) * voronoi_cell)
            i -= 1
        return t - t_new
        
def Lattice(basis=None, inner_product_matrix=None, quadratic_form=None, base_ring=ZZ, **kwargs):
    """
    The `Lattice` function creates lattices using a given base
    or an underlying quadratic form.
    
    INPUT:
    
    - ``basis``
    - ``inner_product_matrix``
    - ``quadratic_form``
    - ``base_ringbase_ring``
    
    OUTPUT:
    
    A lattice.
    
    EXAMPLES::
    
        sage: Lattice([[2, 0, 0], [0, 1, 0]])
        ZZ-lattice of degree 3 and rank 2
        Inner product matrix:
        [1 0]
        [0 4]
        Basis matrix:
        [ 0  1  0]
        [-2  0  0]
        
    A lattice can be specified by a quadratic form::
        
        sage: Lattice(quadratic_form=QuadraticForm(ZZ, 3, [1,2,3,4,5,6]))
        Lattice of degree 3 and rank 3 over Integer Ring
        Inner product matrix:
        [  1   1 3/2]
        [  1   4 5/2]
        [3/2 5/2   6]
        Basis matrix:
        [                  1                   0                   0]
        [                  1  1.732050807568878?                   0]
        [                3/2 0.5773502691896258?  1.848422751068237?]
    """
    if quadratic_form is not None:
        inner_product_matrix = quadratic_form.Gram_matrix_rational()
    if basis is not None:
        if not basis:
            raise ValueError("basis must not be empty")
        #if inner_product_matrix is not None:
        #    raise TypeError("basis and inner_product_matrix cannot both be given")
        basis = matrix(basis)
        rank = basis.dimensions()[0]
        if inner_product_matrix is None:
            inner_product_matrix = basis * basis.transpose()
    elif inner_product_matrix is not None:
        inner_product_matrix = matrix(inner_product_matrix)
        #K = inner_product_matrix.base_ring()
        #rank = len(inner_product_matrix)
        rank = inner_product_matrix.dimensions()[0]
        basis = inner_product_matrix.cholesky()
    else:
        raise TypeError("basis or inner_product_matrix must be given")
    #ambient = base_ring ** rank
    if base_ring == ZZ and inner_product_matrix.base_ring() == ZZ:
        return Lattice_ZZ(base_ring, rank, inner_product_matrix, basis, **kwargs)
    else:
        return Lattice_with_basis(base_ring, rank, inner_product_matrix, basis)
    '''if basis is not None:
        if not basis:
            raise ValueError("basis must not be empty")
        degree = len(basis[0])
        basis_matrix = matrix(basis)
        K = basis_matrix.base_ring()
        if coefficient_ring == ZZ:
            if K <= RR:
                return Lattice_ZZ_in_RR(basis, **kwargs)
            else:
                return Lattice_ZZ(basis)
        else:
            return Lattice_with_basis(coefficient_ring ** degree, basis)
    elif quadratic_form is not None:
        quadratic_form = matrix(quadratic_form)
        # apply Cholesky decomposition to Gram matrix to get basis
        try:
            Q = quadratic_form.cholesky()
        except ValueError:
            # switch to SymbolicRing matrix to allow Cholesky decomposition
            # containing square roots etc.
            dim = quadratic_form.dimensions()
            space = MatrixSpace(SymbolicRing(), dim[0], dim[1])
            quadratic_form = space(quadratic_form)
            Q = quadratic_form.cholesky_decomposition()
        return Lattice(basis=Q, coefficient_ring=coefficient_ring, **kwargs)
    else:
        raise TypeError("basis or quadratic_form must be given")'''
