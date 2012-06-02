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

import sage.matrix.matrix_space

from sage.misc.latex import latex

from sage.modules.free_module import FreeModule_submodule_with_basis_pid
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR

class Lattice_generic(FreeModule_submodule_with_basis_pid):
    """
    Construct a lattice over a PID with a given basis.
    
    INPUT:
    
    - ``ambient`` -- ambient free module over a principal ideal domain `R`,
      i.e. `R^n`;

    - ``basis`` -- list of elements of `K^n`, where `K` is the fraction field
      of `R`. These elements must be linearly independent and will be used as
      the default basis of the constructed submodule.
    """
        
    def __init__(self, ambient, gens):
        """
        See :class:`Lattice_generic` for documentation.
        """
        super(Lattice_generic, self).__init__(ambient, gens, echelonize=False,
            already_echelonized=True, check=False)
        # Don't perform check to keep original base field
        # (e.g. RR^3-vectors in ZZ-lattice, without coercing them to QQ^3).
        # Consequently, the basis has to be coerced to appropriate types before
        
        self._element_class = gens[0].parent()
        
    def basis_matrix(self):
        try:
            return self.__basis_matrix
        except AttributeError:
            # Use self.element_class().base_ring() instead of self.base_ring()
            MAT = sage.matrix.matrix_space.MatrixSpace(self.element_class().base_ring(),
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
    
class Lattice_ZZ(Lattice_generic):
    def __init__(self, gens):
        gens = list(gens)
        degree = len(gens[0])
        super(Lattice_ZZ, self).__init__(ZZ ** degree, gens)
        
class Lattice_ZZ_in_RR(Lattice_ZZ):
    def __init__(self, gens):
        super(Lattice_ZZ_in_RR, self).__init__(gens)
        
def Lattice(basis=None, coefficient_ring=ZZ, quadratic_form=None):
    """
    The `Lattice` function creates lattices using a given base
    or an underlying quadratic form.
    
    INPUT:
    
    -  ``basis``
    -  ``coefficient_ring``
    -  ``quadratic_form``
    
    OUTPUT:
    
    A lattice.
    
    EXAMPLES::
    
        sage: Lattice([[1, 0], [0,1]])
        Lattice of degree 3 and rank 2 over Integer Ring
        Basis matrix:
        [1 0]
        [0 1]
    """
    if basis is not None:
        if not basis:
            raise ValueError("basis must not be empty")
        basis_matrix = matrix(basis)
        degree = len(basis[0])
        K = basis_matrix.base_ring()
        Kn = K ** degree
        basis = [Kn(element) for element in basis]
        if coefficient_ring == ZZ:
            if K == RR:
                return Lattice_ZZ_in_RR(basis)
            else:
                return Lattice_ZZ(basis)
        else:
            return Lattice_generic(coefficient_ring ** degree, basis)
    elif quadratic_form is not None:
        raise NotImplementedError()
    else:
        raise TypeError("base or quadratic_form must be given")
