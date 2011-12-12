"""BoxClaw State class."""

import pyclaw.state


class State(pyclaw.state.State):
    """See the corresponding PyClaw class documentation."""

    # def meqn():
    #     doc = r"""(int) - Number of unknowns (components of q)"""
    #     def fget(self):
    #         if self.q_mf is None:
    #             raise Exception('state.meqn has not been set.')
    #         else: return self.q_mf.dof
    #     return locals()
    # meqn = property(**meqn())

    # def mp():
    #     doc = r"""(int) - Number of derived quantities (components of p)"""
    #     def fset(self,mp):
    #         if self._p_mf is not None:
    #             raise Exception('You cannot change state.mp after p is initialized.')
    #         else:
    #             self._p_mf = self._create_MF(mp)
    #             self.gpVec = self._p_mf.createGlobalVector()
    #     def fget(self):
    #         if self._p_mf is None:
    #             raise Exception('state.mp has not been set.')
    #         else: return self._p_mf.dof
    #     return locals()
    # mp = property(**mp())

    # def mF():
    #     doc = r"""(int) - Number of derived quantities (components of p)"""
    #     def fset(self,mF):
    #         if self._F_mf is not None:
    #             raise Exception('You cannot change state.mp after p is initialized.')
    #         else:
    #             self._F_mf = self._create_MF(mF)
    #             self.gFVec = self._F_mf.createGlobalVector()
    #     def fget(self):
    #         if self._F_mf is None:
    #             raise Exception('state.mF has not been set.')
    #         else: return self._F_mf.dof
    #     return locals()
    # mF = property(**mF())

    # def maux():
    #     doc = r"""(int) - Number of auxiliary fields"""
    #     def fget(self):
    #         if self.aux_mf is None: return 0
    #         else: return self.aux_mf.dof
    #     return locals()
    # maux = property(**maux())

    def q():
        """Array to store solution (q) values.

        Settting state.meqn automatically allocates space for q, as does
        setting q itself.

        Assume we're using the first box of the multifab for now.
        """
        def fget(self):
            if self.q_mf is None: return 0
            fab = self.q_mf.fab(1)
            return fab.array
        def fset(self, q):
            meqn = q.shape[0]
            if self.q_mf is None: self._init_q_mf(meqn)
            fab = self.q_mf.fab(1)
            fab[...] = q
        return locals()
    q = property(**q())

    # def p():
    #     r"""
    #     Array containing values of derived quantities for output.
    #     """
    #     def fget(self):
    #         if self._p_mf is None: return 0
    #         shape = self.grid.ng
    #         shape.insert(0,self.mp)
    #         p=self.gpVec.getArray().reshape(shape, order = 'F')
    #         return p
    #     def fset(self,p):
    #         mp = p.shape[0]
    #         if self.gpVec is None: self.init_p_mf(mp)
    #         self.gpVec.setArray(p.reshape([-1], order = 'F'))
    #     return locals()
    # p = property(**p())

    # def F():
    #     r"""
    #     Array containing pointwise values (densities) of output functionals.
    #     This is just used as temporary workspace before summing.
    #     """
    #     def fget(self):
    #         if self._F_mf is None: return 0
    #         shape = self.grid.ng
    #         shape.insert(0,self.mF)
    #         F=self.gFVec.getArray().reshape(shape, order = 'F')
    #         return F
    #     def fset(self,F):
    #         mF = F.shape[0]
    #         if self.gFVec is None: self.init_F_mf(mF)
    #         self.gFVec.setArray(F.reshape([-1], order = 'F'))
    #     return locals()
    # F = property(**F())

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
            fab = self.aux_mf.fab(1)
            return fab.array
        def fset(self,aux):
            # It would be nice to make this work also for parallel
            # loading from a file.
            if self.aux_mf is None:
                maux = aux.shape[0]
                self._init_aux_mf(maux)
            fab = self.aux_mf.fab(1)
            fab[...] = aux
        return locals()
    aux = property(**aux())

    def ndim():
        def fget(self):
            return self.grid.ndim
        return locals()
    ndim = property(**ndim())


    def __init__(self, grid, meqn, maux=0):
        """

        Here we don't call super because q and aux must be properties in PetClaw
        but should not be properties in PyClaw.

        """

        import boxclaw.grid
        if not isinstance(grid, boxclaw.grid.Grid):
            raise Exception("""A BoxClaw State object must be initialized with
                             a BoxClaw Grid object.""")
        self.aux_mf = None
        self.q_mf = None
        self.p_mf = None
        self.F_mf = None

        self.grid = grid
        self.aux_global = {}
        self.t = 0.0
        self.mcapa = -1

        self._init_q_mf(meqn)

        if maux > 0:
            self._init_aux_mf(maux)

    def _init_aux_mf(self, maux, mbc=0):
        """Initialize BoxLib multifab for the auxiliary array *aux*."""

        self.aux_mf = self._create_multifab(maux, mbc)

        # self.aux_mf = self._create_MF(maux,mbc)
        # self.gauxVec = self.aux_mf.createGlobalVector()
        # self.lauxVec = self.aux_mf.createLocalVector()

    def _init_q_mf(self,meqn,mbc=0):
        """Initialize BoxLib multifab for the solution *q*."""

        self.q_mf = self._create_multifab(meqn, mbc)

        # self.q_mf = self._create_MF(meqn,mbc)
        # self.gqVec = self.q_mf.createGlobalVector()
        # self.lqVec = self.q_mf.createLocalVector()

        # #Now set the local indices for the Dimension objects:
        # ranges = self.q_mf.getRanges()
        # for i,nrange in enumerate(ranges):
        #     dim = self.grid.dimensions[i]
        #     dim.nstart = nrange[0]
        #     dim.nend   = nrange[1]
        #     dim.lowerg = dim.lower + dim.nstart*dim.d

    def _create_multifab(self, nc, mbc=0):
        """Create a BoxLib multifab."""

        from pyboxlib import layout, multifab

        if self.ndim == 1:
            boxes = [ [(1,1,1), (self.grid.n[0],1,1)] ]
        elif self.ndim == 2:
            raise NotImplemented
        elif self.ndim == 3:
            raise NotImplemented
        else:
            raise Exception("Invalid number of dimensions")

        la = layout()
        la.create(boxes=boxes)

        mf = multifab()
        mf.create(la, components=nc, ghost_cells=mbc, interleave=True)

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
