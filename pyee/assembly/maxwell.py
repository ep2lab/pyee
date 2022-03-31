import numpy as np
import scipy.sparse as scp
from . import ampere
from . import ampere2
from . import ampere3
from . import faraday
from . import boundaries

def constructProblem(data, mode):
    # Allocate
    nz = data['geometry']['nz']
    nx = data['geometry']['nx']
    size = 6 * (nz + 1) * (nx + 1)
    row    = []
    column = []
    entry  = []
    b = np.zeros([size, 1], dtype=complex) 


    # Internal Boundaries Flags
    if not('BCs' in data):
        EqEz = np.ones((nz,nx),np.int8)
        EqEx = np.ones((nz,nx),np.int8)
        EqEy = np.ones((nz,nx),np.int8)
    else:
        EqEz = data['BCs']['EqEz']
        EqEx = data['BCs']['EqEx']
        EqEy = data['BCs']['EqEy']

    # Assembly
    row, column, entry, b = ampere3.constructProblem(data, mode, row, column, entry, b, EqEz, EqEx, EqEy)
    row, column, entry, b = faraday.constructProblem(data, mode, row, column, entry, b)
    row, column, entry, b = boundaries.constructProblem(data, mode, row, column, entry, b, EqEz, EqEx, EqEy)

    A = scp.csr_matrix((entry, (row, column)),(size, size), dtype=complex)

    return A, b




