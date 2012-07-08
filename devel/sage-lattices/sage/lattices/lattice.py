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
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.libs.pari.gen import pari
        
from diamond_cutting import calculate_voronoi_cell

class Lattice_with_basis(FreeModule_submodule_with_basis_pid):
    """
    Construct a general lattice over a PID with a given basis.
    
    INPUT:
    
    - ``ambient`` -- ambient free module over a principal ideal domain `R`,
      i.e. `R^n`;

    - ``basis`` -- list of elements of `K^n`, where `K` is the fraction field
      of `R`.
      
    EXAMPLES::
    
        sage: Lattice([[2, 0], [0, 1]], GF(3))
        Lattice of degree 2 and rank 2 over Finite Field of size 3
        Basis matrix:
        [2 0]
        [0 1]
    """
    
    def __init__(self, ambient, basis):
        """
        See :class:`Lattice_with_basis` for documentation.
        
        TESTS::
        
            sage: L = Lattice([[1, 0], [0, 1]])
            sage: L.random_element() # random
            (-1, 2)
            sage: TestSuite(L).run()
        """
        
        super(Lattice_with_basis, self).__init__(ambient, basis, echelonize=False,
            already_echelonized=True)
        
    def _repr_(self):
        """
        Text representation of this lattice.
        """
        return "Lattice of degree %s and rank %s over %s\nBasis matrix:\n%s"%(
            self.degree(), self.rank(), self.base_ring(), self.basis_matrix())
    
class Lattice_ZZ(Lattice_with_basis):
    """
    Construct a lattice with integer coefficients.
    
    EXAMPLES::
    
        sage: Lattice([[GF(3)(1), 0], [0, 1]])
        Integer lattice of degree 2 and rank 2
        Basis matrix:
        [1 0]
        [0 1]
    """
    def __init__(self, basis):
        """
        See :class:`Lattice_ZZ` for documentation.
        
        TESTS::
        
            sage: L = Lattice([[GF(3)(1), 0], [0, 1]])
            sage: TestSuite(L).run()
        """
        basis = list(basis)
        degree = len(basis[0])
        super(Lattice_ZZ, self).__init__(ZZ ** degree, basis)
        
    def _repr_(self):
        """
        Text representation of this lattice.
        """
        return "Integer lattice of degree %s and rank %s\nBasis matrix:\n%s"%(
            self.degree(), self.rank(), self.basis_matrix())
        
class Lattice_ZZ_in_RR(Lattice_ZZ):
    """
    Construct a lattice that is embedded in a real-valued space
    (e.g. RR^n, QQ^n, or ZZ^n).
    
    EXAMPLES::
        
        sage: L = Lattice([[1.0, 0, 0], [0, 1, 0]]); L
        Real-embedded integer lattice of degree 3 and rank 2
        Basis matrix:
        [1 0 0]
        [0 1 0]
        
    By default, the basis is reduced using the LLL algorithm::
    
        sage: L = Lattice([[6, 1], [9, 0]]); L
        Real-embedded integer lattice of degree 2 and rank 2
        Basis matrix:
        [ 0  3]
        [ 3 -1]
        sage: L.discriminant()
        81
        
    However, you can prevent this::
    
        sage: Lattice([[6, 1], [9, 0]], reduce_basis=False)
        Real-embedded integer lattice of degree 2 and rank 2
        Basis matrix:
        [6 1]
        [9 0]
        sage: L.discriminant()
        81
    """
    def __init__(self, basis, reduce_basis=True):
        """
        See :class:`Lattice_ZZ_in_RR` for documentation.
        
        TESTS::
        
            sage: L = Lattice([[2, 1, 0], [0, 1, 0]])
            sage: TestSuite(L).run()
        """
        if reduce_basis:
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
            self.__reduced_basis = None
        super(Lattice_ZZ_in_RR, self).__init__(basis)
        
    def reduced_basis(self):
        if self.__reduced_basis is None:
            basis_matrix = self.basis_matrix()
            basis_matrix = basis_matrix.LLL()
            basis = list(v for v in basis_matrix if v)
            self.__reduced_basis = basis
        return self.__reduced_basis
        
    def _repr_(self):
        """
        Text representation of this lattice.
        """
        return "Real-embedded integer lattice of degree %s and rank %s\nBasis matrix:\n%s"%(
            self.degree(), self.rank(),  self.basis_matrix())
        
    def shortest_vectors(self, max_length=None, max_count=None):
        """
        Find shortest vectors using Pari's Fincke-Pohst algorithm.
        
        sage: L = Lattice([[2, 0], [0, 3]])
        sage: L.shortest_vectors()
        [(2, 0), (0, 3)]
        sage: L = Lattice([[2, 1], [1, 1]])
        sage: L.shortest_vectors()
        [(0, 1), (-1, 0)]
        """
        qf = self.gram_matrix()
        if max_length is None:
            # choose trivial upport bound for vector length
            max_length = sum(sum(x) ** 2 for x in self.basis())
        if max_count is None:
            max_count = self.degree()
        #if self.base_ring() == QQ:
        #    flag = 2
        #else:
        #    flag = 0
        flag = 2 # allow non-integral entries 
        count, length, vectors = pari(qf).qfminim(max_length, max_count, flag)
        vectors = vectors.python()
        return [self.linear_combination_of_basis(v) for v in vectors.columns()]
    
    def voronoi_cell(self, radius=None):
        """
        Compute the Voronoi cell of a lattice, returning a Polyhedron.
        
        INPUT:
        
        - ``radius``  -- radius of ball containing considered vertices
          (default: automatic determination).
          
        OUTPUT:
        
        The Voronoi cell as a Polyhedron instance.
        
        EXAMPLES::
        
            sage: L = Lattice([[1, 0], [0, 1]])
            sage: V = L.voronoi_cell()
            sage: V.Vrepresentation()
            (A vertex at (1/2, -1/2), A vertex at (1/2, 1/2), A vertex at (-1/2, 1/2), A vertex at (-1/2, -1/2))
            
        Lattices not having full dimension are handled as well:
        
            sage: L = Lattice([[2, 0, 0], [0, 2, 0]])
            sage: V = L.voronoi_cell()
            sage: V.Hrepresentation()
            (An inequality (-1, 0, 0) x + 1 >= 0, An inequality (0, -1, 0) x + 1 >= 0, An inequality (1, 0, 0) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0)
        
        "Over-dimensional" lattices are reduced first:
        
            sage: L = Lattice([[1, 0], [2, 0], [0, 2]])
            sage: L.voronoi_cell().Vrepresentation()
            (A vertex at (1/2, -1), A vertex at (1/2, 1), A vertex at (-1/2, 1), A vertex at (-1/2, -1))
        """
        
        basis_matrix = matrix(self.reduced_basis())
        return calculate_voronoi_cell(basis_matrix, radius=radius)
        
def Lattice(basis=None, coefficient_ring=ZZ, quadratic_form=None, **kwargs):
    """
    The `Lattice` function creates lattices using a given base
    or an underlying quadratic form.
    
    INPUT:
    
    - ``basis``
    - ``coefficient_ring``
    - ``quadratic_form``
    
    OUTPUT:
    
    A lattice.
    
    EXAMPLES::
    
        sage: Lattice([[2, 0, 0], [0, 1, 0]])
        Real-embedded integer lattice of degree 3 and rank 2
        Basis matrix:
        [0 1 0]
        [2 0 0]
    """
    if basis is not None:
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
        raise NotImplementedError()
    else:
        raise TypeError("basis or quadratic_form must be given")
