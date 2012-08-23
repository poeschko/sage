"""
Lattice plot functionality

AUTHORS:

- Jan Poeschko (2012-08-15): initial version
"""

from sage import plot

def subsets(s):
    """
    Generate all subsets of a list of elements.
    
    INPUT:
    
    - ``s`` -- a list of elements.
    
    OUTPUT:
    
    A generator iterating through all subsets of `s`.
    """
    n_elems = len(s)
    n_subsets = 2**len(s)
    for i in range(0,n_subsets):
        sub = []
        for j in range(0,n_elems):
            if (i >> j & 1):
                sub.append(s[j])
        yield sub

def plot_lattice(L, **options):
    """
    Plot a lattice.
    
    The sum of each subset of basis vectors (including their negated counterparts)
    is plotted.
    
    INPUT:
    
    - ``L`` -- a `Lattice` instance;
    
    - ``options`` - options to be passed to `list_plot`.
    
    OUTPUT:
    
    Graphics containing the plot of the lattice.
    
    EXAMPLES::
    
        sage: L = Lattice([[2, 0], [0, 3]])
        sage: plot_lattice(L)
        
    Three-dimensional lattices are supported, too::
    
        sage: L = Lattice([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
        sage: plot_lattice(L)
        
    Higher dimensions are not supported, though::
    
        sage: L = random_lattice(4)
        sage: plot_lattice(L)
        Traceback (most recent call last):
        ...
        ValueError: only 2-dimensional and 3-dimensional lattices can be plotted
    """
    dim = L.dimension()
    basis = L.embedded_basis()
    #points = sum(([-v, v] for v in basis), []) + [L(0)]
    #basis = sum(([-v, v] for v in basis), [])
    points = [sum(v, L(0)) for v in subsets(basis)]
    points = [tuple(v.list()) for v in points]
    points = set(points)
    #print sorted(points)
    #return points
    return plot.point.points(points)

    if dim == 2:
        return list_plot(points, **options)
    elif dim == 3:
        return list_plot3d(points, **options)
    else:
        raise ValueError("only 2-dimensional and 3-dimensional lattices can be plotted")
