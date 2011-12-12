"""BoxClaw controller class."""

from pyclaw.controller import Controller as pyclawController

class Controller(pyclawController):
    def __init__(self):
        super(Controller,self).__init__()

        self.output_format = 'ascii'

    def is_proc_0(self):
        from pyboxlib import fboxlib
        rank = fboxlib.mpi_rank()
        return rank == 0
