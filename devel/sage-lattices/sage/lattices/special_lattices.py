"""
Special named lattices

This module provides a constructor for some special named lattices.

AUTHORS:

- Jan Poeschko (2012-08-10): initial version

"""

def special_lattice(name):
    """
    Construct a special named lattice.
    
    INPUT:
    
    - ``name`` -- the name of the lattice to construct (as a string).
      The list of supported names is returned by :func:`special_lattice_names`.
    
    OUTPUT:
    
    A lattice.
    
    EXAMPLES::
    
    Consider the Leech lattice::
        
        sage: leech = special_lattice("Leech")
        sage: leech
        ZZ-lattice of degree 24 and rank 24
        Inner product matrix:
        24 x 24 dense matrix over Integer Ring
        Basis matrix:
        24 x 24 dense matrix over Algebraic Field
        
    The Leech lattice is unimodular::
    
        sage: leech.is_unimodular()
        True
    """
    from constructor import Lattice
    
    basis = _basis.get(name)
    if basis is not None:
        return Lattice(basis=basis, reduce=False)
    ipm = _inner_product_matrices.get(name)
    if ipm is not None:
        return Lattice(inner_product_matrix=ipm)
    raise ValueError("unkown special lattice: %s" % name)

def special_lattice_names():
    """
    Return the list of supported special named lattices.
    
    OUTPUT:
    
    A list of strings containing the supported names.
    
    EXAMPLES::
    
        sage: special_lattice_names()
        ["BodyCenteredCubic"]
    """
    return sorted(_inner_product_matrices.keys() | _basis.keys())

_basis = {
    'BodyCenteredCubic': [
        [2, 0, 0],
        [0, 2, 0],
        [1, 1, 1]
    ],
    'FaceCenteredCubic': [
        [-1, -1, 0],
        [1, -1, 0],
        [0, 1, -1]
    ],
    'HexagonalLattice': [
        [-1, 1],
        [0, -1]
    ],
    'SimpleCubic': [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ],
    'SquareLattice': [
        [1, 0],
        [0, 1]
    ],
}

_inner_product_matrices = {
    'Leech': [
        [8, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 2, 4, 2, 2, 2, 4, 2, 2, 2, 0, 0, 0, -3],
        [4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 1, 2, 1, 0, 0, -1],
        [4, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 2, 2, 1, 1, 1, 0, 0, -1],
        [4, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 0, 0, -1],
        [4, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, 0, -1],
        [4, 2, 2, 2, 2, 4, 2, 2, 2, 2, 2, 1, 2, 2, 1, 1, 2, 1, 2, 1, 0, 0, 0, -1],
        [4, 2, 2, 2, 2, 2, 4, 2, 2, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 0, 0, 0, -1],
        [2, 2, 2, 2, 2, 2, 2, 4, 1, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 0, 0, 1],
        [4, 2, 2, 2, 2, 2, 2, 1, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, -1],
        [4, 2, 2, 2, 2, 2, 2, 1, 2, 4, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 0, 1, 0, -1],
        [4, 2, 2, 2, 2, 2, 2, 1, 2, 2, 4, 2, 2, 1, 2, 1, 2, 1, 2, 1, 0, 0, 1, -1],
        [2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 4, 1, 2, 2, 2, 1, 2, 2, 2, 2, 1, 1, 1],
        [4, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 4, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, -1],
        [2, 2, 1, 1, 2, 2, 1, 2, 2, 2, 1, 2, 2, 4, 2, 2, 1, 2, 2, 2, 2, 2, 1, 1],
        [2, 1, 2, 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 4, 2, 1, 2, 2, 2, 2, 1, 2, 1],
        [2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 4, 1, 2, 2, 2, 2, 1, 1, 1],
        [4, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 4, 2, 2, 2, 1, 1, 1, -1],
        [2, 1, 2, 1, 2, 1, 1, 2, 2, 2, 1, 2, 1, 2, 2, 2, 2, 4, 2, 2, 2, 2, 1, 1],
        [2, 1, 1, 2, 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 4, 2, 2, 1, 2, 1],
        [2, 2, 1, 1, 2, 1, 2, 2, 2, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 4, 2, 1, 1, 1],
        [0, 1, 1, 1, 1, 0, 0, 2, 1, 0, 0, 2, 1, 2, 2, 2, 1, 2, 2, 2, 4, 2, 2, 2],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 4, 2, 2],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 2, 4, 2],
        [-3, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 2, 2, 2, 4]
    ],
    'SimpleHexagonal': [
        [1, 0, 0],
        [0, 2, 1],
        [0, 1, 2]
    ],
    'SimpleOrthorhombic': [
        [1, 0, 0],
        [0, 2, 0],
        [0, 0, 3]
    ],
}
