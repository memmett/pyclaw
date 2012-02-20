
import numpy as np
from ctypes import *

so = CDLL('./libsharpclaw_.so')

class WKSPACE(Structure):
    _fields_ = [
        ('num_eqn', c_int),
        ('num_dim', c_int),
        ('num_waves', c_int), 
        ('index_capa', c_int),
        ('maxmx', c_int),
        ('num_ghost', c_int),
        ('char_decomp', c_int),
        ('lim_type', c_int),
        ('multid_recon', c_int), 
        ('weno_order', c_int),
        ('fwave', c_int),
        ('tfluct_solver', c_int),
        ('epweno', c_double),
        ('bptr', c_void_p) ]


class PROBLEMDATA(Structure):
    _fields_ = [
        ('name', c_char_p),
        ('data', POINTER(c_double)) ] 

class RP(Structure):
    _fields_ = [
        ('context', c_void_p),
        ('num_pdata', c_int),
        ('pdata', POINTER(PROBLEMDATA)*32) ]


class wrapper(object):
    """Python wrapper to the libsharpclaw1.so shared library.

    >>> sc = sharpclaw1()
    >>> sc.params.maxnx = 320

    """

    def __init__(self):
        """Create a new workspace."""

        self.rp = None
        self.pfloat64 = POINTER(c_double)

        self.cptr = c_void_p()
        so.create_workspace(byref(self.cptr))

        self.params = WKSPACE.from_address(self.cptr.value)



    def setup(self):
        """Allocate work space arrays etc."""

        so.setup_workspace(byref(self.cptr))


    def set_dx(self, dx):
        """Set dx."""

        dx = np.array(dx, order='F')
        so.set_dx(self.cptr, dx.ctypes.data_as(self.pfloat64), c_int(dx.shape[0]))


    def set_problem_data(self, problem_data):

        rp = RP()
        rp.num_pdata = 0

        so.rp_register(byref(rp))

        for i in range(rp.num_pdata):
            key = str(rp.pdata[i].contents.name)
            if key:
                value = problem_data[key]
                if isinstance(value, float):
                    d = c_double(value)
                    rp.pdata[i].contents.data.contents.value = value
                else:
                    raise TypeError("Type of problem data '%s' not supported yet." % key)
            else:
                raise ValueError("Problem data '%s' not found." % key)

        self.rp = rp


    def flux1(self, q, aux, dt, t, mx):
        """Call flux1."""

        dq  = np.zeros((self.params.num_eqn,self.params.maxmx+2*self.params.num_ghost),
                       order='F')
        cfl = c_double()

        so.flux1(self.cptr,
                 q.ctypes.data_as(self.pfloat64), 
                 dq.ctypes.data_as(self.pfloat64),
                 c_double(dt), 
                 byref(cfl),
                 c_double(t), 
                 c_int(1), 
                 self.params.num_eqn,
                 c_int(mx),
                 self.params.num_ghost,
                 self.params.maxmx,
                 byref(self.rp))

        return dq, cfl.value

    def teardown(self):
        """XXX"""
        so.rp_options_teardown(byref(self.rp))
        so.rp_teardown(byref(self.rp))

        self.rp = None


if __name__ == '__main__':
    import doctest
    doctest.testmod()
