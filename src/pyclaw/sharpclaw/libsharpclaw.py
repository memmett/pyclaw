
import numpy as np
from ctypes import *

so = CDLL('./libsharpclaw_.so')

class WKSPACE(Structure):
    _fields_ = [
        ('maxmx', c_int),
        ('num_ghost', c_int),
        ('num_eqn', c_int),
        ('num_dim', c_int),
        ('num_waves', c_int), 
        ('index_capa', c_int),
        ('char_decomp', c_int),
        ('lim_type', c_int),
        ('multid_recon', c_int), 
        ('weno_order', c_int),
        ('fwave', c_int),
        ('tfluct_solver', c_int),
        ('epweno', c_double) ]


class wrapper(object):
    """Python wrapper to the libsharpclaw1.so shared library.

    >>> sc = sharpclaw1()
    >>> sc.params.maxnx = 320

    """

    def __init__(self):
        """Create a new workspace."""

        self.cptr = c_void_p()
        so.create_workspace(byref(self.cptr))
        self.params = WKSPACE.from_address(self.cptr.value)

        self.pfloat64 = POINTER(c_double)

    def setup(self):
        """Allocate work space arrays etc."""

        so.setup_workspace(byref(self.cptr))

    def flux1(self, q, aux, dt, t, mx):
        """Call flux1."""

        dq  = np.empty(q.shape)
        cfl = np.empty(q.shape[0])

        print 'flux1...'
        so.flux1(self.cptr,
                 q.ctypes.data_as(self.pfloat64), 
                 dq.ctypes.data_as(self.pfloat64),
                 aux.ctypes.data_as(self.pfloat64), 
                 c_float(dt), 
                 cfl.ctypes.data_as(self.pfloat64),
                 c_float(t), 
                 c_int(1), 
                 c_int(mx),
                 self.params.num_ghost,
                 self.params.maxmx)
        print 'flux1... done.'

        print dq

        return dq, cfl


if __name__ == '__main__':
    import doctest
    doctest.testmod()
