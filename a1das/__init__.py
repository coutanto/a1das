#
# Various python modules to read/reduct/process Febus A1 DAS file and
# Febus A1 DAS socket stream.
#
# core.py : basic IO operation on files and socket stream. Definition of basic Classes
#     _core_febus_file.py : for Febus file
#     _core_file_reducted.py : for reducted_filefiles
#     _core_socket.py : for socket
#
# reduction.py : reduction and conversion operations from strain or raw data
#     Reduction to simple hdf5 structures with/without transposition, decimation
#     antialiasing filter.
#     Conversion from raw optical phase to strain performed by a binary
#     module converted from fortran by f2py
#     raw2strPy.???.so : binary for conversion
#
# xcor.py: computation of cross-correlations from strain file or socket stream
#     Manage the configuration of the cross-correlation and the processing parameters
#
# a1das version 2.0.0

__all__=["core", "reduction", "xcor", "util", "_plot.py", "_a1das_exception.py", "_core_reducted_file.py",
         "_a1xcorPy", "_sosfilterPy", "_core_febus_file.py", "vserver"]
#
# from a1das import *
# importera core, reduction, xcor, plot
#
#
# import a1das
# importera toutes les declarations de classes, procedures de core
# visible sans mettre le module core
#
# ET imortera tous les modules reduction, xcor, ...
from .core import *

from . import reduction
from . import xcor
from . import util
from . import vserver
from . import _core_febus_file
from . import _core_reducted_file
from . import _core_socket
from . import _a1xcorPy
from . import _sosfilterPy
from . import _plot
