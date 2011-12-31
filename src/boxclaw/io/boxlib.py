"""Routines for reading and writing BoxLib-style output files."""

import os,sys
import logging

import pyclaw.solution

logger = logging.getLogger('io')

import pickle
import numpy as np

from pyboxlib import fboxlib, layout, multifab


def write_boxlib(solution, frame,
                 path='./', file_prefix='claw',
                 write_aux=False, write_p=False,
                 options={}):
    """Write pickle and BoxLib data files representing the solution."""

    defaults = { 'clobber': True }
    defaults.update(options)
    options = defaults

    clobber = options['clobber']

    pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    boxlib_filename = os.path.join(path, '%s.bxl' % file_prefix) + str(frame).zfill(4)

    if solution.maux > 0 and write_aux:
        write_aux = True
        aux_filename = os.path.join(path, '%s_aux.bxl' % file_prefix)
    else:
        write_aux = False

    if not clobber:
        if os.path.exists(pickle_filename):
            raise IOError('Cowardly refusing to clobber %s!' % pickle_filename)
        if os.path.exists(boxlib_filename):
            raise IOError('Cowardly refusing to clobber %s!' % boxlib_filename)
        if write_aux and os.path.exists(aux_filename):
            raise IOError('Cowardly refusing to clobber %s!' % aux_filename)

    rank = fboxlib.mpi_rank()

    if rank == 0:
        pickle_file = open(pickle_filename,'wb')
        if write_p:
            pickle.dump({'t':solution.t,'meqn':solution.mp,'nstates':len(solution.states),
                         'maux':solution.maux,'ndim':solution.ndim,'write_aux':write_aux,
                         'aux_global' : solution.aux_global}, pickle_file)
        else:
            pickle.dump({'t':solution.t,'meqn':solution.meqn,'nstates':len(solution.states),
                         'maux':solution.maux,'ndim':solution.ndim,'write_aux':write_aux,
                         'aux_global' : solution.aux_global}, pickle_file)

    # write
    for state in solution.states:
        grid = state.grid

        if rank == 0:
            pickle.dump({'level':grid.level,
                         'names':grid.name,'lower':grid.lower,
                         'n':grid.n,'d':grid.d}, pickle_file)

        if write_p:
            state.p_mf.write(boxlib_filename, 'p')
        else:
            state.q_mf.write(boxlib_filename, 'q')

        if write_aux:
            state.aux_mf.write(aux_filename, 'aux')

    if rank == 0:
        pickle_file.close()

# def read_petsc(solution,frame,path='./',file_prefix='claw',read_aux=False,options={}):
#     r"""
#     Read in pickles and PETSc data files representing the solution

#     :Input:
#      - *solution* - (:class:`~pyclaw.solution.Solution`) Solution object to
#        read the data into.
#      - *frame* - (int) Frame number to be read in
#      - *path* - (string) Path to the current directory of the file
#      - *file_prefix* - (string) Prefix of the files to be read in.
#        ``default = 'fort'``
#      - *read_aux* (bool) Whether or not an auxiliary file will try to be read
#        in.  ``default = False``
#      - *options* - (dict) Optional argument dictionary, see
#        `PETScIO Option Table`_

#     .. _`PETScIO Option Table`:

#     format   : one of 'ascii' or 'binary'

#     """

#     # Option parsing
#     option_defaults = {'format':'binary'}

#     for (k,v) in option_defaults.iteritems():
#         if options.has_key(k):
#             pass
#         else:
#             options[k] = option_defaults[k]

#     pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
#     viewer_filename = os.path.join(path, '%s.ptc' % file_prefix) + str(frame).zfill(4)
#     aux_viewer_filename = os.path.join(path, '%s_aux.ptc' % file_prefix)
#     if frame < 0:
#         # Don't construct file names with negative frameno values.
#         raise IOError("Frame " + str(frame) + " does not exist ***")

#     pickle_file = open(pickle_filename,'rb')

#     # this dictionary is mostly holding debugging information, only nstates is needed
#     # most of this information is explicitly saved in the individual grids
#     value_dict = pickle.load(pickle_file)
#     nstates    = value_dict['nstates']
#     ndim       = value_dict['ndim']
#     maux       = value_dict['maux']
#     meqn       = value_dict['meqn']

#     # now set up the PETSc viewer
#     if options['format'] == 'ascii':
#         viewer = PETSc.Viewer().createASCII(viewer_filename, PETSc.Viewer.Mode.READ)
#         if read_aux:
#             aux_viewer = PETSc.Viewer().createASCII(aux_viewer_filename, PETSc.Viewer.Mode.READ)
#     elif options['format'] == 'binary':
#         if hasattr(PETSc.Viewer,'createMPIIO'):
#             viewer = PETSc.Viewer().createMPIIO(viewer_filename, PETSc.Viewer.Mode.READ)
#         else:
#             viewer = PETSc.Viewer().createBinary(viewer_filename, PETSc.Viewer.Mode.READ)
#         if read_aux:
#             if os.path.exists(aux_viewer_filename):
#                 if hasattr(PETSc.Viewer,'createMPIIO'):
#                     aux_viewer = PETSc.Viewer().createMPIIO(aux_viewer_filename, PETSc.Viewer.Mode.READ)
#                 else:
#                     aux_viewer = PETSc.Viewer().createBinary(aux_viewer_filename, PETSc.Viewer.Mode.READ)
#             else:
#                 from warnings import warn
#                 aux_file_path = os.path.join(path,aux_viewer_filename)
#                 warn('read_aux=True but aux file %s does not exist' % aux_file_path)
#                 read_aux=False
#     else:
#         raise IOError('format type %s not supported' % options['format'])

#     for m in xrange(nstates):
#         grid_dict = pickle.load(pickle_file)

#         level   = grid_dict['level']
#         names   = grid_dict['names']
#         lower   = grid_dict['lower']
#         n       = grid_dict['n']
#         d       = grid_dict['d']

#         import petclaw ##
#         dimensions = []
#         for i in xrange(ndim):
#             dimensions.append(
#                 #pyclaw.solution.Dimension(names[i],lower[i],lower[i] + n[i]*d[i],n[i]))
#                 petclaw.Dimension(names[i],lower[i],lower[i] + n[i]*d[i],n[i]))
#         #grid = pyclaw.solution.Grid(dimensions)
#         grid = petclaw.Grid(dimensions)
#         grid.level = level
#         #state = pyclaw.state.State(grid)
#         state = petclaw.State(grid,meqn,maux) ##
#         state.t = value_dict['t']
#         state.aux_global = value_dict['aux_global']

# #       DA View/Load is broken in Petsc-3.1.8, we can load/view the DA if needed in petsc-3.2
# #       state.q_da.load(viewer)
#         state.gqVec.load(viewer)

#         if read_aux:
#             state.gauxVec.load(aux_viewer)

#         solution.states.append(state)

#     pickle_file.close()
#     viewer.destroy()
#     if read_aux:
#         aux_viewer.destroy()

# def read_petsc_t(frame,path='./',file_prefix='claw'):
#     r"""Read only the petsc.pkl file and return the data

#     :Input:
#      - *frame* - (int) Frame number to be read in
#      - *path* - (string) Path to the current directory of the file
#      - *file_prefix* - (string) Prefix of the files to be read in.
#        ``default = 'claw'``

#     :Output:
#      - (list) List of output variables
#       - *t* - (int) Time of frame
#       - *meqn* - (int) Number of equations in the frame
#       - *ngrids* - (int) Number of grids
#       - *maux* - (int) Auxillary value in the frame
#       - *ndim* - (int) Number of dimensions in q and aux

#     """

#     base_path = os.path.join(path,)
#     path = os.path.join(base_path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
#     try:
#         f = open(path,'rb')
#         logger.debug("Opening %s file." % path)
#         grid_dict = pickle.load(f)

#         t      = grid_dict['t']
#         meqn   = grid_dict['meqn']
#         nstates = grid_dict['nstates']
#         maux   = grid_dict['maux']
#         ndim   = grid_dict['ndim']

#         f.close()
#     except(IOError):
#         raise
#     except:
#         logger.error("File " + t_fname + " should contain t, meqn, ngrids, maux, ndim")
#         print "File " + t_fname + " should contain t, meqn, ngrids, maux, ndim"
#         raise

#     return t,meqn,nstates,maux,ndim
