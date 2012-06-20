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
from sage.rings.real_mpfr import RR

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
    """
    def _init_(self, basis):
        """
        See :class:`Lattice_ZZ_in_RR` for documentation.
        
        TESTS::
        
            sage: L = Lattice([[2, 1, 0], [0, 1, 0]])
            sage: TestSuite(L).run()
        """
        super(Lattice_ZZ_in_RR, self).__init__(basis)
        
    def _repr_(self):
        """
        Text representation of this lattice.
        """
        return "Real-embedded integer lattice of degree %s and rank %s\nBasis matrix:\n%s"%(
            self.degree(), self.rank(),  self.basis_matrix())
        
def Lattice(basis=None, coefficient_ring=ZZ, quadratic_form=None):
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
        [2 0 0]
        [0 1 0]
    """
    if basis is not None:
        if not basis:
            raise ValueError("basis must not be empty")
        basis_matrix = matrix(basis)
        degree = len(basis[0])
        K = basis_matrix.base_ring()
        if coefficient_ring == ZZ:
            if K <= RR:
                return Lattice_ZZ_in_RR(basis)
            else:
                return Lattice_ZZ(basis)
        else:
            return Lattice_with_basis(coefficient_ring ** degree, basis)
    elif quadratic_form is not None:
        raise NotImplementedError()
    else:
        raise TypeError("basis or quadratic_form must be given")
