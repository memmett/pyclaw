#!/usr/bin/env python
# encoding: utf-8
"""
1D acoustics example.
"""

import os, sys

try:
    import numpy as np
    from petsc4py import PETSc
    
except:
    sys.path.append("/opt/share/ksl/petsc4py/dev-aug29/ppc450d/lib/python/")
    sys.path.append("/opt/share/ksl/numpy/dev-aug29/ppc450d/lib/python/")
    
    import numpy as np
    from petsc4py import PETSc

from petclaw.grid import PCGrid as Grid
from petclaw.grid import PCDimension as Dimension
from pyclaw.solution import Solution
from petclaw.evolve.petclaw import PetClawSolver1D
from pyclaw.controller import Controller

def qinit(grid):
    grid.init_q_petsc_structures()

    # Initial Data parameters
    ic = 1
    beta = 100.
    gamma = 0.
    x0 = 0.75
    x1 = 0.7
    x2 = 0.9
    
    x =grid.x.center
    q=np.zeros([len(x),grid.meqn], order = 'F')
    
    # Gaussian
    qg = np.exp(-beta * (x-x0)**2) * np.cos(gamma * (x - x0))
    # Step Function
    qs = (x > x1) * 1.0 - (x > x2) * 1.0
    
    if ic == 1: q[:,0] = qg
    elif ic == 2: q[:,0] = qs
    elif ic == 3: q[:,0] = qg + qs
    q[:,1]=0.
    grid.q=q


# Initialize grids and solutions
x = Dimension('x',0.0,1.0,100,mthbc_lower=1,mthbc_upper=1)
grid = Grid(x)
grid.aux_global['zz']=1.0
grid.aux_global['cc']=1.0
grid.meqn=2
grid.t = 0.0
qinit(grid)
init_solution = Solution(grid)

# Solver setup
solver = PetClawSolver1D(kernelsType = 'F')

solver.dt = 0.0004
solver.max_steps = 5000
solver.set_riemann_solver('acoustics')
solver.order = 2
solver.mthlim = [4,4]

useController = True
makePlot = True


if useController:

    # Controller instantiation
    claw = Controller()
    claw.outdir = './_output/'
    claw.keep_copy = True
    claw.nout = 50
    claw.outstyle = 1
    claw.output_format = 'petsc'
    claw.tfinal = 1.0
    claw.solutions['n'] = init_solution
    claw.solver = solver

    # Solve
    status = claw.run()

    if makePlot:
        if claw.keep_copy:
    
            for n in xrange(0,11):
                
                sol = claw.frames[n]
                plotTitle="time: {0}".format(sol.t)
                viewer = PETSc.Viewer.DRAW(sol.grid.gqVec.comm)
                OptDB = PETSc.Options()
                OptDB['draw_pause'] = 1
                viewer(sol.grid.gqVec)

else:
    sol = {"n":init_solution}
    solver.evolve_to_time(sol,0.25)

    sol = sol["n"]

    if makePlot:
        viewer = PETSc.Viewer.DRAW(grid.gqVec.comm)
        OptDB = PETSc.Options()
        OptDB['draw_pause'] = -1
        viewer(grid.gqVec)
