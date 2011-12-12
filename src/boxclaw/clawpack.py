"""BoxClaw solvers."""

import pyclaw.clawpack

class ClawSolver1D(pyclaw.clawpack.ClawSolver1D):
    """BoxClaw solver for 1D problems using classic Clawpack algorithms."""
    pass

class ClawSolver2D(pyclaw.clawpack.ClawSolver2D):
    """BoxClaw solver for 2D problems using classic Clawpack algorithms."""
    pass

class ClawSolver3D(pyclaw.clawpack.ClawSolver3D):
    """PetClaw solver for 3D problems using classic Clawpack algorithms."""
    pass
