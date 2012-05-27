"""
Lattices

This module provides the abstract base class :class:`Lattice_generic` from which
all lattices in Sage derive, as well as a selection of more
specific base classes.

The class inheritance hierarchy is:

- :class:`Lattice_generic`

  - :class:`LatticeZZ`

    - :class:`Lattice_embedded_in_RR`

AUTHORS:

- Jan Poeschko (2012-05-26): initial version.

"""

from sage.modules.free_module import FreeModule_ambient

class Lattice_generic(FreeModule_ambient):
    def __init__(self, base_ring, base):
        FreeModule_ambient.__init__(self, base_ring, len(base))
        self.base = base
        
    def _repr_(self):
        """
        Text representation of this lattice.
        """
        return "Lattice on %s with base %s" % (self.base_ring(), self.base)
    
    def _latex_(self):
        """
        LaTeX representation of this lattice.
        """
        return "\\text{Lattice}"
    
    def _matrix_(self):
        """
        Matrix representation of this lattice.
        """
        return self.base
        
def Lattice(base=None, quadratic_form=None):
    """
    The `Lattice` function creates lattices using a given base
    or an underlying quadratic form.
    
    INPUT:
    
    -  ``base``
    -  ``quadratic_form``
    
    OUTPUT:
    
    A lattice.
    
    EXAMPLES::
    
        sage: Lattice([1, 2])
        Lattice on Integer Ring with base [1, 2]
    """
    
    if base is not None:
        if base:
            return Lattice_generic(base[0].parent(), base)
        else:
            raise ValueError("base must not be empty")
    elif quadratic_form is not None:
        raise NotImplementedError()
    else:
        raise TypeError("base or quadratic_form must be given")
