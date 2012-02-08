#!/usr/bin/env python
# encoding: utf-8

import numpy as np

nx = 40

def mapc2p(patch, x):

    print 'x', x

    dx = 1.0/nx

    m = x
    m[5] -= 0.1*dx
    m[6] += 0.1*dx
    m[7] -= 0.1*dx
    m[8] += 0.1*dx
    m[9] -= 0.1*dx
    m[10] += 0.1*dx

    # print 'm', m

    return m

def acoustics(use_petsc=False,
              kernel_language='Fortran',
              solver_type='classic',
              weno_order=5,
              mapped=False,
              unifcoeffs=False,
              iplot=False,
              htmlplot=False,
              outdir='./_output'):
    """
    This example solves the 1-dimensional acoustics equations in a homogeneous
    medium.
    """
    from numpy import sqrt, exp, cos

    #=================================================================
    # Import the appropriate classes, depending on the options passed
    #=================================================================
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver1D()
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
        solver.weno_order=weno_order
    else: raise Exception('Unrecognized value of solver_type.')

    #========================================================================
    # Instantiate the solver and define the system of equations to be solved
    #========================================================================
    solver.kernel_language=kernel_language
    from clawpack.riemann import rp_acoustics
    solver.num_waves=rp_acoustics.num_waves

    if kernel_language=='Python': 
        solver.rp = rp_acoustics.rp_acoustics_1d
    else:
        from clawpack.riemann import rp1_acoustics
        solver.rp = rp1_acoustics

    solver.limiters = pyclaw.limiters.tvd.MC

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.aux_bc_lower[0] = pyclaw.BC.periodic
    solver.aux_bc_upper[0] = pyclaw.BC.periodic

    #========================================================================
    # Instantiate the domain and set the boundary conditions
    #========================================================================
    x  = pyclaw.Dimension('x', 0.0, 1.0, nx)
    domain = pyclaw.Domain(x)
    num_eqn = 2
    state = pyclaw.State(domain,num_eqn)

    if mapped:
        dx = 1.0/nx
        num_ghost = (weno_order+1)/2
        # state.grid.mapc2p = lambda p, x: x + 0.8*dx*np.sin(2*np.pi*x)
        state.grid.mapc2p = mapc2p
        state.grid.compute_p_edges(recompute=True)
        state.grid.compute_p_centers(recompute=True)
        if not unifcoeffs:
            print 'precomputing weno coeffs...'
            state.precompute_mapped_weno(solver.weno_order)
            print 'precomputing weno coeffs... done.'
        area = state.grid.p_edges[0][1:] - state.grid.p_edges[0][:-1]
        print area
        state.aux = np.zeros((1,nx))
        state.aux[0,:] = area/dx
        state.index_capa = 0

    #========================================================================
    # Set problem-specific variables
    #========================================================================
    rho = 1.0
    bulk = 1.0

    state.problem_data['rho']=rho
    state.problem_data['bulk']=bulk
    state.problem_data['zz']=sqrt(rho*bulk)
    state.problem_data['cc']=sqrt(bulk/rho)
 

    #========================================================================
    # Set the initial condition
    #========================================================================
    xc=domain.grid.x.centers
    beta=100; gamma=0; x0=0.75
    if mapped:
        from scipy.integrate import quadrature
        q0 = lambda xc: exp(-beta * (xc-x0)**2) * cos(gamma * (xc - x0))
        x = state.grid.p_edges[0]
        for i in range(state.q.shape[1]):
            a = x[i]
            b = x[i+1]
            print a, b, b-a
            q, e = quadrature(q0, a, b)
            state.q[0,i] = q/(b-a)
    else:
        state.q[0,:] = exp(-beta * (xc-x0)**2) * cos(gamma * (xc - x0))
    q1 = np.array(state.q[0,:])

    state.q[1,:] = 0.

    solver.dt_initial=domain.grid.delta[0]/state.problem_data['cc']*0.1

    #========================================================================
    # Set up the controller object
    #========================================================================
    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.keep_copy = True
    claw.num_output_times = 5

    claw.tfinal = 1.0

    # Solve
    status = claw.run()

    q2 = state.q[0,:]
    print abs(q1-q2).max()

    # Plot results
    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics)
