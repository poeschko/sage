"""
all.py -- much of sage is imported into this module, so you don't
          have to import everything individually.
    
TESTS:
    
    This is to test #10570. If the number of stackframes at startup
    changes due to a patch you made, please check that this was an
    intended effect of your patch.

    ::
        
        sage: import gc
        sage: import inspect
        sage: from sage import *
        sage: frames=[x for x in gc.get_objects() if inspect.isframe(x)]
        sage: len(frames)
        11
    
"""

###############################################################################
#
#   SAGE: System for Algebra and Geometry Experimentation    
#
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
###############################################################################

# Error message that matches the Sage/IPython defaults
quit = "Use Ctrl-D (i.e. EOF), %Exit, or %Quit to exit without confirmation."
exit = quit

import os, sys

if 'SAGE_ROOT' not in os.environ:
    raise RuntimeError("To use the Sage libraries, set the environment variable SAGE_ROOT to the Sage build directory and LD_LIBRARY_PATH to $SAGE_ROOT/local/lib")
 
if sys.version_info[:2] < (2, 5):
    print >>sys.stderr, "Sage requires Python 2.5 or newer"
    sys.exit(1)

###################################################################

# We have to set this here so urllib, etc. can detect it.
import sage.server.notebook.gnutls_socket_ssl
sage.server.notebook.gnutls_socket_ssl.require_SSL()

###################################################################

from sage.ext.c_lib import _init_csage, sig_on_count
_init_csage()

from time                import sleep

from sage.misc.all       import *         # takes a while

from sage.misc.sh import sh

from sage.libs.all       import *

from sage.rings.all      import *
from sage.matrix.all     import *

# This must come before Calculus -- it initializes the Pynac library.
import sage.symbolic.pynac

from sage.modules.all    import *
from sage.monoids.all    import *
from sage.algebras.all   import *
from sage.modular.all    import *
from sage.schemes.all    import *
from sage.graphs.all     import *
from sage.groups.all     import *
from sage.databases.all  import *
from sage.structure.all  import *
from sage.categories.all import *
from sage.sets.all       import *
from sage.probability.all import *
from sage.interfaces.all import *

from sage.symbolic.all   import *

from sage.functions.all  import *
from sage.calculus.all   import *

from sage.server.all     import *
import sage.tests.all as tests

from sage.crypto.all     import *
import sage.crypto.mq as mq

from sage.plot.all       import *
from sage.plot.plot3d.all     import *

from sage.coding.all     import *
from sage.combinat.all   import *

from sage.lfunctions.all import *

from sage.geometry.all   import *
from sage.geometry.triangulation.all   import *

from sage.homology.all   import *

from sage.quadratic_forms.all import *

from sage.gsl.all        import *

from sage.games.all      import *

from sage.media.all      import *

from sage.logic.all      import *

from sage.numerical.all  import *

from sage.stats.all      import *
import sage.stats.all as stats

import sage.finance.all  as finance

import sage.interacts.all as interacts
from sage.interacts.debugger import debug 

from sage.parallel.all   import *

from sage.ext.fast_callable  import fast_callable
from sage.ext.fast_eval      import fast_float

from sage.sandpiles.all import *

from sage.tensor.all     import *

from sage.lattices.all import *

from copy import copy, deepcopy

# The code executed here uses a large amount of Sage components
from sage.rings.qqbar import _init_qqbar
_init_qqbar()

#Deprecate the is_* functions from the top level
#All of these functions should be removed from the top level
#after a few releases, and this code should be removed.
#--Mike Hansen 9/25/2008
message = "\nUsing %(name)s from the top level is deprecated since it was designed to be used by developers rather than end users.\nIt most likely does not do what you would expect it to do.  If you really need to use it, import it from the module that it is defined in."
sage.misc.superseded.deprecated_callable_import(
    10107, None, globals(), locals(),
    [name for name in globals().keys() if name.startswith('is_') and name[3].isupper()],
    message)

del message, name


###########################################################
#### WARNING:
# DO *not* import numpy / matplotlib / networkx here!!
# Each takes a surprisingly long time to initialize,
# and that initialization should be done more on-the-fly
# when they are first needed.
###########################################################

###################################################################

# maximize memory resources
#try:
#    import resource   # unix only...
#    resource.setrlimit(resource.RLIMIT_AS, (-1,-1))
#except:
#    pass

# very useful 2-letter shortcuts
CC = ComplexField()
QQ = RationalField()
RR = RealField()  # default real field
ZZ = IntegerRing()
# NOTE: QQ, RR, and ZZ are used by the pre-parser, and should not be
# overwritten by the user, unless they want to change the meaning of
# int and real in the interpreter (which is a potentially valid thing
# to do, and doesn't mess up anything else in the Sage library).
# E.g., typing "int = ZZ" in the Sage interpreter makes int literals
# acts as Python ints again.



# Some shorter shortcuts:
# Q = QQ
# Z = ZZ
# C = CC
#i = CC.gen(0)
true = True
false = False

oo = infinity
#x = PolynomialRing(QQ,'x').gen()

from sage.misc.copying import license
copying = license
copyright = license

_cpu_time_ = cputime()
_wall_time_ = walltime()

def quit_sage(verbose=True):
    """
    If you use Sage in library mode, you should call this function
    when your application quits.

    It makes sure any child processes are also killed, etc.
    """
    if verbose:
        t1 = cputime(_cpu_time_)
        t1m = int(t1/60); t1s=t1-t1m*60
        t2 = walltime(_wall_time_)
        t2m = int(t2/60); t2s=t2-t2m*60
        print "Exiting Sage (CPU time %sm%.2fs, Wall time %sm%.2fs)."%(
               t1m,t1s,t2m,t2s)
    from sage.interfaces.quit import expect_quitall
    expect_quitall(verbose=verbose)

    import sage.matrix.matrix_mod2_dense
    sage.matrix.matrix_mod2_dense.free_m4ri()

    import sage.libs.flint.flint
    sage.libs.flint.flint.free_flint_stack()

    pari._unsafe_deallocate_pari_stack()
    
    ### The following is removed -- since it would cleanup
    ### the tmp directory that the sage cleaner depends upon.
    # The following code close all open file descriptors,
    # so that on shared file systems the delete_tmpfiles
    # command below works.
    # AUTHOR:
    #    * Kate Minola (2007-05-03)
    #import resource             # Resource usage information.
    #maxfd = resource.getrlimit(resource.RLIMIT_NOFILE)[1]
    #if maxfd != resource.RLIM_INFINITY:
        # Iterate through and close all file descriptors.
    #    for fd in range(0, maxfd):
    #        try:
    #            os.close(fd)
    #        except OSError:  # ERROR, fd wasn't open to begin with (ignored)
    #            pass
    # Now delete the temp files
    #from sage.misc.misc import delete_tmpfiles
    #delete_tmpfiles()

    # stop the twisted reactor
    try:
       from twisted.internet import reactor
       if reactor.running:
          reactor.callFromThread(reactor.stop)
    except ImportError:
       pass

    # Free globally allocated mpir integers.
    import sage.rings.integer
    sage.rings.integer.free_integer_pool()
    sage.rings.integer.clear_mpz_globals()
    import sage.algebras.quatalg.quaternion_algebra_element
    sage.algebras.quatalg.quaternion_algebra_element._clear_globals()

    from sage.libs.all import symmetrica 
    symmetrica.end() 
        
def _quit_sage_(self):
    import sage.misc.preparser_ipython
    if sage.misc.preparser_ipython.interface != None:
        sage.misc.preparser_ipython.switch_interface('sage')
        self.exit_now = False
        return
    
    from IPython.genutils import ask_yes_no
    if self.rc.confirm_exit:
        if ask_yes_no('Do you really want to exit ([y]/n)?','y'):
            self.exit_now = True
    else:
        self.exit_now = True
    if self.exit_now:
        quit_sage()
        self.exit_now = True

    return self.exit_now

from IPython.iplib import InteractiveShell
InteractiveShell.exit = _quit_sage_

import sage.misc.displayhook
sage.misc.displayhook.install()

from sage.ext.interactive_constructors_c import inject_on, inject_off

sage.structure.sage_object.register_unpickle_override('sage.categories.category', 'Sets', Sets)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'HeckeModules', HeckeModules)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'Objects', Objects)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'Rings', Rings)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'Fields', Fields)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'VectorSpaces', VectorSpaces)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'Schemes_over_base', sage.categories.schemes.Schemes_over_base)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'ModularAbelianVarieties', ModularAbelianVarieties)
#sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', '', )

# Cache the contents of star imports.
import sage.misc.lazy_import
sage.misc.lazy_import.save_cache_file()


### Debugging for Singular, see trac #10903
# from sage.libs.singular.ring import poison_currRing
# sys.settrace(poison_currRing)


# Write a file indicating that Sage was started up successfully.
def _write_started_file():
    """
    Write a ``sage-started.txt`` file if it does not exist.  The
    contents of this file do not matter, only its existence.
    
    The non-existence of this file will be used as a trigger to run
    ``sage-starts`` during the Sage build.

    TESTS:

    Check that the file exists when Sage is running::

        sage: started_file = os.path.join(SAGE_ROOT, 'local', 'lib', 'sage-started.txt')
        sage: os.path.isfile(started_file)
        True
    """
    started_file = os.path.join(SAGE_ROOT, 'local', 'lib', 'sage-started.txt')
    # Do nothing if the file already exists
    if os.path.isfile(started_file):
        return

    # Current time with a resolution of 1 second
    import datetime
    t = datetime.datetime.now().replace(microsecond=0)

    O = open(started_file, 'w')
    O.write("Sage %s was started at %s\n"%(sage.version.version, t))
    O.close()

_write_started_file()


# Set a new random number seed as the very last thing
# (so that printing initial_seed() and using that seed
# in set_random_seed() will result in the same sequence you got at
# Sage startup).
set_random_seed()
