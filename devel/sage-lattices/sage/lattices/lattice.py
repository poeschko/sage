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

from sage.misc.latex import latex

from sage.modules.free_module import FreeModule_submodule_with_basis_pid
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR

class Lattice_generic(FreeModule_submodule_with_basis_pid):
    def __init__(self, ambient, gens, coefficient_ring):
        super(Lattice_generic, self).__init__(ambient, gens, echelonize=False)
        self._coefficient_ring = coefficient_ring
        
    def _repr_(self):
        """
        Text representation of this lattice.
        """
        return "%s-Lattice in %s" % (self._coefficient_ring, self.vector_space())
    
    def _latex_(self):
        """
        LaTeX representation of this lattice.
        """
        return "\\left\\{\\sum_i a_i v_i\\mid a_i\\in %s,v_i\\in %s\\right\\}" % (
            latex(self._coefficient_ring), 
            '\\left\\{%s\\right\\}' % ','.join(latex(gen) for gen in self.gens()))
    
    def _matrix_(self):
        """
        Matrix representation of this lattice.
        """
        return matrix(self.gens())
    
class Lattice_ZZ(Lattice_generic):
    def __init__(self, ambient, gens):
        super(Lattice_ZZ, self).__init__(ambient, gens, ZZ)
        
class Lattice_ZZ_in_RR(Lattice_ZZ):
    def __init__(self, gens):
        dim = len(gens[0])
        super(Lattice_ZZ_in_RR, self).__init__(RR ** dim, gens)
        
def Lattice(vector_space=None, basis=None, coefficient_ring=ZZ, quadratic_form=None):
    """
    The `Lattice` function creates lattices using a given base
    or an underlying quadratic form.
    
    INPUT:
    
    -  ``vector_space``
    -  ``basis``
    -  ``coefficient_ring``
    -  ``quadratic_form``
    
    OUTPUT:
    
    A lattice.
    
    EXAMPLES::
    
        sage: Lattice([1, 2])
        Lattice on Integer Ring with base [1, 2]
    """
    if basis is not None:
        if not basis:
            raise ValueError("base must not be empty")
        if vector_space is None:
            basis_matrix = matrix(basis)
            vector_space = basis_matrix.base_ring() ** basis_matrix.dimensions()[1]
        if coefficient_ring == ZZ:
            if vector_space.base_ring() == RR:
                return Lattice_ZZ_in_RR(basis)
            else:
                return Lattice_ZZ(vector_space, basis)
        else:
            return Lattice_generic(vector_space, basis, coefficient_ring)
    elif quadratic_form is not None:
        raise NotImplementedError()
    else:
        raise TypeError("base or quadratic_form must be given")
