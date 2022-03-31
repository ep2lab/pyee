import numpy as np
import scipy.sparse as scp

def constructProblem(data, mode, row, column, entry, b, EqEz, EqEx, EqEy):

    """ Constructs the matrix and vector of the linear system """
    axi = data['geometry']['axi']
    freq = data['simulation']['freq'] * 1e6
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

    def kappa(iz, ix):
        kappa = data['kappa'][int(2*iz), int(2*ix), :, :]
        # kappa = data['kappa'][2*int(iz), 2*int(ix), :, :]
        return kappa

    def ja(iz, ix):
        ja = data['current'][int(2*iz), int(2*ix), :]
        return ja

    # ------------------------------ INNER DOMAIN EQUATIONS ------------------------------
    for ix in range(1, nx):
        for iz in range(1, nz):
            ib = diz * iz + dix * ix
        
            # 0 -> Ez
            if EqEz[iz,ix]:
                row.extend([ib + 0] * 10)
                column.append(ib             + 0); entry.append(1j * r(ix) * kappa(iz + 0.5, ix)[0, 0] * dx)
                column.append(ib             + 1); entry.append(1j * r(ix + 0.5) * kappa(iz, ix + 0.5)[0, 1] * dx / 4)
                column.append(ib + diz       + 1); entry.append(1j * r(ix + 0.5) * kappa(iz + 1, ix + 0.5)[0, 1] * dx / 4)
                column.append(ib       - dix + 1); entry.append(1j * r(ix - 0.5) * kappa(iz, ix - 0.5)[0, 1] * dx / 4)
                column.append(ib + diz - dix + 1); entry.append(1j * r(ix - 0.5) * kappa(iz + 1, ix - 0.5)[0, 1] * dx / 4)
                column.append(ib             + 2); entry.append(1j * r(ix) * kappa(iz, ix)[0, 2] * dx / 2)
                column.append(ib + diz       + 2); entry.append(1j * r(ix) * kappa(iz + 1, ix)[0, 2] * dx / 2)
                column.append(ib             + 4); entry.append(-1j * mode * dx)
                column.append(ib             + 5); entry.append(r(ix + 0.5))
                column.append(ib - dix       + 5); entry.append(-r(ix - 0.5))
                b[ib + 0] = ja(iz + 0.5, ix)[0] * r(ix) * dx
            else:
                row.append(ib + 0); column.append(ib + 0); entry.append(1)
                b[ib + 0] = 0

            # 1 -> Ex
            if EqEx[iz,ix]:
                row.extend([ib + 1] * 10)
                column.append(ib             + 0); entry.append(1j * r(ix) * kappa(iz + 0.5, ix)[1, 0] * dz / 4)
                column.append(ib - diz       + 0); entry.append(1j * r(ix) * kappa(iz - 0.5, ix)[1, 0] * dz / 4)
                column.append(ib       + dix + 0); entry.append(1j * r(ix + 1) * kappa(iz + 0.5, ix + 1)[1, 0] * dz / 4)
                column.append(ib - diz + dix + 0); entry.append(1j * r(ix + 1) * kappa(iz - 0.5, ix + 1)[1, 0] * dz / 4)
                column.append(ib             + 1); entry.append(1j * r(ix + 0.5) * kappa(iz, ix + 0.5)[1, 1] * dz)
                column.append(ib             + 2); entry.append(1j * r(ix) * kappa(iz, ix)[1, 2] * dz / 2)
                column.append(ib       + dix + 2); entry.append(1j * r(ix + 1) * kappa(iz, ix + 1)[1, 2] * dz / 2)
                column.append(ib             + 3); entry.append(1j * mode * dz)
                column.append(ib             + 5); entry.append(-r(ix + 0.5))
                column.append(ib - diz       + 5); entry.append(r(ix + 0.5))
                b[ib + 1] = ja(iz, ix + 0.5)[1] * r(ix + 0.5) * dz
            else:
                row.append(ib + 1); column.append(ib + 1); entry.append(1)
                b[ib + 1] = 0

            # 2 -> Ey
            if EqEy[iz,ix]:
                row.extend([ib + 2] * 9)
                column.append(ib             + 0); entry.append(1j * kappa(iz + 0.5, ix)[2, 0] * dz * dx / 2)
                column.append(ib - diz       + 0); entry.append(1j * kappa(iz - 0.5, ix)[2, 0] * dz * dx / 2)
                column.append(ib             + 1); entry.append(1j * kappa(iz, ix + 0.5)[2, 1] * dz * dx / 2)
                column.append(ib       - dix + 1); entry.append(1j * kappa(iz, ix - 0.5)[2, 1] * dz * dx / 2)
                column.append(ib             + 2); entry.append(1j * kappa(iz, ix)[2, 2] * dz * dx)
                column.append(ib             + 3); entry.append(-dz)
                column.append(ib       - dix + 3); entry.append(dz)
                column.append(ib             + 4); entry.append(dx)
                column.append(ib - diz       + 4); entry.append(-dx)
                b[ib + 2] = ja(iz, ix)[2] * dz * dx
            else:
                row.append(ib + 2); column.append(ib + 2); entry.append(1)
                b[ib + 2] = 0


    return row, column, entry, b