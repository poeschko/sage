"""
Lattices

This module provides the abstract base class :class:`Lattice_generic` from which
all lattices in Sage derive, as well as a selection of more
specific base classes.

The class inheritance hierarchy is:

- :class:`Lattice_generic`

  - :class:`Lattice_ZZ`

    - :class:`Lattice_ZZ_in_RR`

AUTHORS:

- Jan Poeschko (2012-05-26): initial version.

"""

#from sage.rings.integer_ring import ZZ

import sage.matrix.matrix_space

from sage.misc.latex import latex

from sage.modules.free_module import FreeModule_submodule_with_basis_pid, element_class
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR

class Lattice_with_basis(FreeModule_submodule_with_basis_pid):
    """
    Construct a lattice over a PID with a given basis.
    
    INPUT:
    
    - ``ambient`` -- ambient free module over a principal ideal domain `R`,
      i.e. `R^n`;

    - ``basis`` -- list of elements of `K^n`, where `K` is the fraction field
      of `R`. These elements must be linearly independent and will be used as
      the default basis of the constructed submodule.
    """
        
    def __init__(self, ambient, basis):
        """
        See :class:`Lattice_with_basis` for documentation.
        """
        #self._Lattice_element_class = basis[0].parent()
        self._Lattice_K = basis[0].parent().base_ring()
        self._Lattice_element_class = element_class(self._Lattice_K, is_sparse=False)
        
        #self.__ambient_module = ambient
        #C = self.element_class()
        #basis = [C(self, x.list(), coerce=False, copy=True) for x in basis]
        #[C(self, x.list(), coerce=False, copy=True) for x in basis]
        
        super(Lattice_with_basis, self).__init__(ambient, basis, echelonize=False,
            already_echelonized=True, check=False)
        # Don't perform check to keep original base field
        # (e.g. RR^3-vectors in ZZ-lattice, without coercing them to QQ^3).
        # Consequently, the basis has to be coerced to appropriate types before
        
    def element_class(self):
        return self._Lattice_element_class
        
    def basis_matrix(self):
        try:
            return self.__basis_matrix
        except AttributeError:
            # Use self._Lattice_K instead of self.base_ring()
            MAT = sage.matrix.matrix_space.MatrixSpace(self._Lattice_K,
                            len(self.basis()), self.degree(),
                            sparse = self.is_sparse())
            if self.is_ambient():
                A = MAT.identity_matrix()
            else:
                A = MAT(self.basis())
            A.set_immutable()
            self.__basis_matrix = A
            return A
        
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
    
        sage: Lattice([[i, 0], [0, 1]])
        Integer lattice of degree 2 and rank 2
        Basis matrix:
        [I 0]
        [0 1]
    """
    def __init__(self, basis):
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
    """
    def __init__(self, basis):
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
        #Kn = K ** degree
        #C = element_class(K, is_sparse=False)
        V = K ** degree
        basis = [V(element) for element in basis]
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
