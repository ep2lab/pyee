import numpy as np

def generate(data):

    zz = 50
    xx = 50
    z_max = 1
    r_max = 1
    mode  = 0 

    E1 = np.zeros([zz, xx, 3], dtype=complex)

    R, Z = np.meshgrid(np.linspace(0, z_max, zz), np.linspace(0, r_max, xx))


    size = np.shape(R)
        
    for i in range(0, size[0]):
        for j in range(0, size[1]):
            rho = R[i, j]
            z   = Z[i, j]
            E1[i, j, 0] = np.sin(np.pi*(1 - rho))*np.sin(np.pi*(1 - z))
            E1[i, j, 1] = np.sin(np.pi*(1 - rho))*np.sin(np.pi*(1 - z))
            E1[i, j, 2] = np.sin(np.pi*(1 - rho))*np.sin(np.pi*(1 - z))

    # dEzz = np.diff(E1[:, :, 0], axis = 0)[:, 1 : zz]
    dErz = np.diff(E1[:, :, 1], axis = 1)[1 : zz, :]
    dEpz = np.diff(E1[:, :, 2], axis = 1)[1 : zz, :]

    dEzr = np.diff(E1[:, :, 0], axis = 0)[:, 1 : xx]
    # dErr = np.diff(E1[:, :, 1], axis = 1)[1 : xx, :]
    dEpr = np.diff(R * E1[:, :, 2], axis = 0)[:, 1 : xx]

    E15 = np.zeros([zz - 1, xx - 1, 3], dtype=complex)
    E15[:, :, 0] = E1[1 : zz, 1 : xx, 0]
    E15[:, :, 1] = E1[1 : zz, 1 : xx, 1]
    E15[:, :, 2] = E1[1 : zz, 1 : xx, 2]

    B1 = np.zeros([zz - 1, xx - 1, 3], dtype=complex)

    B1[:, :, 0] = 1/R[1 : zz, 1 : xx] * (dEpr - 1j*mode*E15[:, :, 1])
    B1[: ,:, 1] = 1/R[1 : zz, 1 : xx] * 1j*mode*E15[:, :, 0] - dEpz
    B1[:, :, 2] = dErz - dEzr

    # dBzz = np.diff(B1[:, :, 0], axis = 0)[:, 1 : zz - 1]
    dBrz = np.diff(B1[:, :, 1], axis = 1)[1 : zz - 1, :]
    dBpz = np.diff(B1[:, :, 2], axis = 1)[1 : zz - 1, :]

    dBzr = np.diff(B1[:, :, 0], axis = 0)[:, 1 : xx - 1]
    # dBrr = np.diff(B1[:, :, 1], axis = 1)[1 : xx - 1, :]
    dBpr = np.diff(R[1 : zz, 1 : xx] * B1[:, :, 2], axis = 0)[:, 1 : xx - 1]


    E = np.zeros([zz - 2, xx - 2, 3], dtype=complex)
    D = E
    E[:, :, 0] = E1[1 : zz - 1, 1 : xx - 1, 0]
    E[:, :, 1] = E1[1 : zz - 1, 1 : xx - 1, 1]
    E[:, :, 2] = E1[1 : zz - 1, 1 : xx - 1, 2]

    B = np.zeros([zz - 2, xx - 2, 3], dtype=complex)
    B[:, :, 0] = B1[0 : zz - 2, 0 : xx - 2, 0]
    B[:, :, 1] = B1[0 : zz - 2, 0 : xx - 2, 1]
    B[:, :, 2] = B1[0 : zz - 2, 0 : xx - 2, 2]

    for i in range(0, zz - 2):
        for j in range(0, xx - 2):
            D[i, j, 0] = data['kappa'][i, j, 0, 0]*E[i, j, 0] + data['kappa'][i, j, 0, 1]*E[i, j, 1] + data['kappa'][i, j, 0, 2]*E[i, j, 2]
            D[i, j, 1] = data['kappa'][i, j, 1, 0]*E[i, j, 0] + data['kappa'][i, j, 1, 1]*E[i, j, 1] + data['kappa'][i, j, 1, 2]*E[i, j, 2]
            D[i, j, 2] = data['kappa'][i, j, 2, 0]*E[i, j, 0] + data['kappa'][i, j, 2, 1]*E[i, j, 1] + data['kappa'][i, j, 2, 2]*E[i, j, 2]

    J = np.zeros([zz - 2, xx - 2, 3], dtype=complex)
    J[:, :, 0] = 1/R[1 : zz - 1, 1 : xx - 1] * (dBpr - 1j*mode*B[:, :, 1]) - 1j*D[:, :, 0]
    J[: ,:, 1] = 1/R[1 : zz - 1, 1 : xx - 1] * 1j*mode*B[:, :, 0] - dBpz   - 1j*D[:, :, 1]
    J[:, :, 2] = dBrz - dBzr - 1j*mode*D[:, :, 2]
    
    return E, J, R[1 : zz - 1, 1 : xx - 1], Z[1 : zz - 1, 1 : xx - 1]



def interpolate(data):
    from scipy.interpolate import interp2d
    E, J, R, Z = generate(data)

    Zsol = data['geometry']['Lz'] * Z
    Rsol = data['geometry']['Lx'] * R

    dz = data['geometry']['Lz'] / data['geometry']['nz']
    dx = data['geometry']['Lx'] / data['geometry']['nx']

    R, Z = np.meshgrid(np.linspace(0, data['geometry']['Lz'], 2*data['geometry']['nz'] + 1), 
                       np.linspace(0, data['geometry']['Lx'], 2*data['geometry']['nx'] + 1))

    ja = np.zeros([np.shape(R)[0], np.shape(R)[1], 3], dtype = complex)

    fzr = interp2d(Zsol, Rsol, np.real(J[:, :, 0]))
    fzi = interp2d(Zsol, Rsol, np.imag(J[:, :, 0]))
    frr = interp2d(Zsol, Rsol, np.real(J[:, :, 1]))
    fri = interp2d(Zsol, Rsol, np.imag(J[:, :, 1]))
    fpr = interp2d(Zsol, Rsol, np.real(J[:, :, 2]))
    fpi = interp2d(Zsol, Rsol, np.imag(J[:, :, 2]))
    for iz in range(0, 2*data['geometry']['nz'] + 1):
        for ix in range(0, 2*data['geometry']['nx'] + 1):
            rho = R[iz, ix]
            z   = Z[iz, ix]
            ja[iz, ix, 0] = fzr(z + dz, rho)[0] + 1j*fzi(z + dz, rho)[0]
            ja[iz, ix, 1] = frr(z, rho + dx)[0] + 1j*fri(z, rho + dx)[0]
            ja[iz, ix, 2] = fpr(z, rho)[0] + 1j*fpi(z, rho)[0]
    return ja