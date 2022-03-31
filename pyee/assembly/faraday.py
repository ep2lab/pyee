import numpy as np
import scipy.sparse as scp

def constructProblem(data, mode, row, column, entry, b):

    """ Constructs the matrix and vector of the linear system """
    axi = data['geometry']['axi']
    freq = data['simulation']['freq']*1e6
    Lz = data['geometry']['Lz'] * 2*np.pi*freq/3e8 # Non-dimensional
    Lx = data['geometry']['Lx'] * 2*np.pi*freq/3e8
    nz = data['geometry']['nz']
    nx = data['geometry']['nx']
    dz = Lz / nz
    dx = Lx / nx
    diz = 6
    dix = 6 * (nz + 1)

    def r(ix):
        if axi:
            r = ix * dx
        else:
            r = 1
        return r

    # ------------------------------ FARADAY'S LAW ------------------------------

    # 3 -> Ez
    for ix in range(0, nx):
        for iz in range(0, nz + 1):
            ib = diz * iz + dix * ix  # block index
            row.extend([ib + 3] * 4)
            column.append(ib       + 3); entry.append(1j * r(ix + 0.5) * dx)
            column.append(ib       + 1); entry.append(1j * mode * dx)
            column.append(ib       + 2); entry.append(r(ix))
            column.append(ib + dix + 2); entry.append(-r(ix + 1))
            b[ib + 3] = 0
    ix = nx  # upper phantom boundary
    for iz in range(0, nz + 1):
        ib = diz * iz + dix * ix  # block index
        row.append(ib + 3); column.append(ib + 3); entry.append(1)
        b[ib + 3] = 0

    # 4 -> Ex
    for ix in range(1, nx + 1):
        for iz in range(0, nz):
            ib = diz * iz + dix * ix  # block index
            row.extend([ib + 4] * 4)
            column.append(ib       + 4); entry.append(1j * r(ix) * dz)
            column.append(ib       + 0); entry.append(-1j * mode * dz)
            column.append(ib       + 2); entry.append(-r(ix))
            column.append(ib + diz + 2); entry.append(r(ix))
            b[ib + 4] = 0
    iz = nz  # right phantom boundary
    for ix in range(0, nx + 1):
        ib = diz * iz + dix * ix  # block index
        row.append(ib + 4); column.append(ib + 4); entry.append(1)
        b[ib + 4] = 0
    ix = 0  # lower boundary
    for iz in range(0, nz):
        ib = diz * iz + dix * ix  # block index
        if axi == 1 and (mode == 1 or mode == -1):
            row.extend([ib + 4] * 4)
            column.append(ib       + 4); entry.append(-1j * dz * dx)
            column.append(ib + dix + 0); entry.append(1j * mode * dz)
            column.append(ib       + 2); entry.append(dx)
            column.append(ib + diz + 2); entry.append(- dx)
            b[ib + 4] = 0
        else:  # PEC or mode != 1,-1
            row.append(ib + 4); column.append(ib + 4); entry.append(1)
            b[ib + 4] = 0
    # 5 -> Ey
    for ix in range(0, nx):
        for iz in range(0, nz):
            ib = diz * iz + dix * ix  # block index
            row.extend([ib + 5] * 5)
            column.append(ib             + 5); entry.append(1j * dz * dx)
            column.append(ib             + 0); entry.append(-dz)
            column.append(ib       + dix + 0); entry.append(dz)
            column.append(ib             + 1); entry.append(dx)
            column.append(ib + diz       + 1); entry.append(-dx)
            b[ib + 5] = 0
    ix = nx  # upper phantom boundary
    for iz in range(0, nz + 1):
        ib = diz * iz + dix * ix  # block index
        row.append(ib + 5); column.append(ib + 5); entry.append(1)
        b[ib + 5] = 0
    iz = nz  # right phantom boundary
    for ix in range(0, nx):
        ib = diz * iz + dix * ix  # block index
        row.append(ib + 5); column.append(ib + 5); entry.append(1)
        b[ib + 5] = 0
    # Return matrix and vector

    return row, column, entry, b