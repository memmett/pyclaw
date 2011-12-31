"""BoxClaw State class.

This class follows the PetClaw state class fairly closely.

"""

import pyclaw.state
import boxclaw.grid

from pyboxlib import fboxlib, layout, multifab


class State(pyclaw.state.State):

    def __init__(self, grid, meqn, maux=0):

        if not isinstance(grid, boxclaw.grid.Grid):
            raise Exception("""A BoxClaw State object must be initialized with
                             a BoxClaw Grid object.""")
        self.aux_mf = None
        self.q_mf = None
        self.p_mf = None
        self.F_mf = None

        self.layout = None
        self.box = None

        self.grid = grid
        self.aux_global = {}
        self.t = 0.0
        self.mcapa = -1

        self._init_q_mf(meqn)

        if maux > 0:
            self._init_aux_mf(maux)


    def ndim():
        def fget(self):
            return self.grid.ndim
        return locals()
    ndim = property(**ndim())


    def meqn():
        doc = r"""(int) - Number of unknowns (components of q)"""
        def fget(self):
            if self.q_mf is None:
                raise Exception('state.meqn has not been set.')
            else: return self.q_mf.nc
        return locals()
    meqn = property(**meqn())


    def mp():
        doc = r"""(int) - Number of derived quantities (components of p)"""
        def fset(self,mp):
            if self.p_mf is not None:
                raise Exception('You cannot change state.mp after p is initialized.')
            else:
                self.p_mf = self._create_multifab(mp)
        def fget(self):
            if self.p_mf is None:
                raise Exception('state.mp has not been set.')
            else: return self.p_mf.nc
        return locals()
    mp = property(**mp())


    def mF():
        doc = r"""(int) - Number of derived quantities (components of F)"""
        def fset(self,mF):
            if self.F_mf is not None:
                raise Exception('You cannot change state.mF after F is initialized.')
            else:
                self.F_mf = self._create_multifab(mF)
        def fget(self):
            if self.F_mf is None:
                raise Exception('state.mF has not been set.')
            else: return self.F_mf.nc
        return locals()
    mF = property(**mF())


    def maux():
        doc = r"""(int) - Number of auxiliary fields"""
        def fget(self):
            if self.aux_mf is None: return 0
            else: return self.aux_mf.nc
        return locals()
    maux = property(**maux())


    def q():
        """Array to store solution (q) values.

        Settting state.meqn automatically allocates space for q, as does
        setting q itself.
        """
        def fget(self):
            if self.q_mf is None: return 0
            fab = self.q_mf.fab(self.box)
            return fab.array
        def fset(self, q):
            meqn = q.shape[0]
            if self.q_mf is None: self._init_q_mf(meqn)
            fab = self.q_mf.fab(self.box)
            fab[...] = q
        return locals()
    q = property(**q())


    def p():
        """Array containing values of derived quantities for output."""
        def fget(self):
            if self.p_mf is None: return 0
            fab = self.p_mf.fab(self.box)
            return fab.array
        def fset(self,p):
            if self.p_mf is None: self._init_p_mf()
            fab = self.p_mf.fab(self.box)
            fab[...] = p
        return locals()
    p = property(**p())


    def F():
        """Array containing pointwise values (densities) of output functionals.

        This is just used as temporary workspace before summing.
        """
        def fget(self):
            if self.F_mf is None: return 0
            fab = self.F_mf.fab(self.box)
            return fab.array
        def fset(self,F):
            mF = F.shape[0]
            if self.F_mf is None: self._init_F_mf()
            fab = self.p_F.fab(self.box)
            fab[...] = F
        return locals()
    F = property(**F())


    def aux():
        """
        We never communicate aux values; every processor should set its own ghost cell
        values for the aux array.  The global aux vector is used only for outputting
        the aux values to file; everywhere else we use the local vector.
        """
        def fget(self):
            if self.aux_mf is None: return None
            shape = self.grid.ng
            shape.insert(0,self.maux)
            fab = self.aux_mf.fab(self.box)
            return fab.array
        def fset(self,aux):
            if self.aux_mf is None:
                maux = aux.shape[0]
                self._init_aux_mf(maux)
            fab = self.aux_mf.fab(self.box)
            fab[...] = aux
        return locals()
    aux = property(**aux())


    def _init_q_mf(self,meqn,mbc=0):
        """Initialize BoxLib multifab for the solution *q*."""

        self.q_mf = self._create_multifab(meqn, mbc)


    def _init_p_mf(self):
        raise NotImplemented


    def _init_F_mf(self):
        raise NotImplemented


    def _init_aux_mf(self, maux, mbc=0):
        """Initialize BoxLib multifab for the auxiliary array *aux*."""

        self.aux_mf = self._create_multifab(maux, mbc)


    def _create_layout(self):

        # create boxes to divide domain amoungst processors
        boxes = []
        if self.ndim == 1:
            n = self.grid.n[0] / fboxlib.mpi_size()
            for k in range(fboxlib.mpi_size()):
                boxes.append( ((k*n+1,1,1), ((k+1)*n,1,1)) )
        elif self.ndim == 2:
            n = self.grid.n[1] / fboxlib.mpi_size()
            for k in range(fboxlib.mpi_size()):
                boxes.append( ((1,k*n+1,1), (1,(k+1)*n,1)) )
        elif self.ndim == 3:
            n = self.grid.n[2] / fboxlib.mpi_size()
            for k in range(fboxlib.mpi_size()):
                boxes.append( ((1,1,k*n+1), (1,1,(k+1)*n)) )
        else:
            raise Exception("Invalid number of dimensions.")

        # create layout
        la = layout()
        la.create(boxes=boxes)
        self.box = la.local_boxes[0]
        self.layout = la

        assert len(la.local_boxes) == 1

        # set local indices of the Dimension objects
        if self.ndim == 1:
            n = self.grid.n[0] / fboxlib.mpi_size()
            k = self.box-1

            dim = self.grid.dimensions[0]
            dim.nstart = k*n+1
            dim.nend   = (k+1)*n+1

        elif self.ndim == 2:
            raise NotImplemented
        elif self.ndim == 3:
            raise NotImplemented
        else:
            raise Exception("Invalid number of dimensions.")


    def _create_multifab(self, nc, mbc=0):
        """Create a BoxLib multifab."""

        if self.layout is None:
            self._create_layout()

        mf = multifab()
        mf.create(self.layout, components=nc, ghost_cells=mbc, interleave=True)

        return mf


    # def set_q_from_qbc(self,mbc,qbc):
    #     """
    #     Set the value of q using the array qbc. for PetSolver, this
    #     involves setting qbc as the local vector array then perform
    #     a local to global communication.
    #     """

    #     grid = self.grid
    #     if grid.ndim == 1:
    #         self.q = qbc[:,mbc:-mbc]
    #     elif grid.ndim == 2:
    #         self.q = qbc[:,mbc:-mbc,mbc:-mbc]
    #     elif grid.ndim == 3:
    #         self.q = qbc[:,mbc:-mbc,mbc:-mbc,mbc:-mbc]
    #     else:
    #         raise NotImplementedError("The case of 3D is not handled in "\
    #         +"this function yet")

    # def get_qbc_from_q(self,mbc,whichvec,qbc):
    #     """
    #     Returns q with ghost cells attached.  For PetSolver,
    #     this means returning the local vector.
    #     """
    #     shape = [n + 2*mbc for n in self.grid.ng]

    #     if whichvec == 'q':
    #         self.q_mf.globalToLocal(self.gqVec, self.lqVec)
    #         shape.insert(0,self.meqn)
    #         return self.lqVec.getArray().reshape(shape, order = 'F')

    #     elif whichvec == 'aux':
    #         self.aux_mf.globalToLocal(self.gauxVec, self.lauxVec)
    #         shape.insert(0,self.maux)
    #         return self.lauxVec.getArray().reshape(shape, order = 'F')

    # def set_mbc(self,mbc):
    #     r"""
    #     This is a hack to deal with the fact that petsc4py
    #     doesn't allow us to change the stencil_width (mbc).

    #     Instead, we initially create DAs with stencil_width=0.
    #     Then, in solver.setup(), we call this function to replace
    #     those DAs with new ones that have the right stencil width.

    #     This could be made more efficient using some PETSc calls,
    #     but it only happens once so it seems not to be worth it.
    #     """
    #     q0 = self.q.copy()
    #     self._init_q_mf(self.meqn,mbc)
    #     self.q = q0

    #     if self.aux is not None:
    #         aux0 = self.aux.copy()
    #         self._init_aux_mf(self.maux,mbc)
    #         self.aux = aux0

    # def sum_F(self,i):
    #     return self.gFVec.strideNorm(i,0)
