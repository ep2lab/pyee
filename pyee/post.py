import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interpn
from numpy import logical_and as AND
import scipy.ndimage.filters as scp
from matplotlib import path
from matplotlib import ticker, cm
import h5py
import os
from matplotlib.colors import LogNorm

def z(data, iz):
    dz = data['geometry']['Lz'] / data['geometry']['nz']
    z = iz * dz
    return z


def x(data, ix):
    dx = data['geometry']['Lx'] / data['geometry']['nx']
    x = ix * dx
    return x


def r(data, ix):
    dx = data['geometry']['Lx'] / data['geometry']['nx']
    if data['geometry']['axi']:
        r = ix * dx
    else:
        r = 1
    return r

def getField(data, solution, field):

    """ From a solution vector returns the mesh and the field """
    nz = data['geometry']['nz']
    nx = data['geometry']['nx']
    vz = np.linspace(0, nz, nz + 1)
    vx = np.linspace(0, nx, nx + 1)
    F = np.zeros([nz + 1, nx + 1], dtype = complex)
    Z = np.zeros([nz + 1, nx + 1])
    R = np.zeros([nz + 1, nx + 1])

    """ Factor for magnetic field """
    FAC   = (3e8)**(-1)

    if field == 'Ez':
        zz = z(data, vz[0:-1] + 0.5)
        xx = x(data, vx)
        R, Z = np.meshgrid(xx, zz)
        F = np.zeros([nz, nx + 1], dtype = complex)
        for ix in range(0, nx + 1):
            for iz in range(0, nz):
                F[iz, ix] = solution[6 * iz + 6 * (nz + 1) * ix + 0]
    elif field == 'Ey':
        zz = z(data, vz)
        xx = x(data, vx)
        R, Z = np.meshgrid(xx, zz)
        F = np.zeros([nz + 1, nx + 1], dtype = complex)
        for ix in range(0, nx + 1):
            for iz in range(0, nz + 1):
                F[iz, ix] = solution[6 * iz + 6 * (nz + 1) * ix + 2]
    elif field == 'Ex':
        zz = z(data, vz)
        xx = x(data, vx[0:-1] + 0.5)
        R, Z = np.meshgrid(xx, zz)
        F = np.zeros([nz + 1, nx], dtype = complex)
        for ix in range(0, nx):
            for iz in range(0, nz + 1):
                F[iz, ix] = solution[6 * iz + 6 * (nz + 1) * ix + 1]
    elif field == 'Bz':
        zz = z(data, vz)
        xx = x(data, vx + 0.5)
        R, Z = np.meshgrid(xx, zz)
        F = np.zeros([nz + 1, nx + 1], dtype = complex)
        for ix in range(0, nx + 1):
            for iz in range(0, nz + 1):
                F[iz, ix] = solution[6 * iz + 6 * (nz + 1) * ix + 3]*FAC
    elif field == 'Bx':
        zz = z(data, vz + 0.5)
        xx = x(data, vx)
        R, Z = np.meshgrid(xx, zz)
        F = np.zeros([nz + 1, nx + 1], dtype = complex)
        for ix in range(0, nx + 1):
            for iz in range(0, nz + 1):
                F[iz, ix] = solution[6 * iz + 6 * (nz + 1) * ix + 4]*FAC
    elif field == 'By':
        zz = z(data, vz + 0.5)
        xx = x(data, vx + 0.5)
        R, Z = np.meshgrid(xx, zz)
        F = np.zeros([nz + 1, nx + 1], dtype = complex)
        for ix in range(0, nx + 1):
            for iz in range(0, nz + 1):
                F[iz, ix] = solution[6 * iz + 6 * (nz + 1) * ix + 5]*FAC
    else:
        pass

    return Z, R, F

def interpFields(data, sol, field_list = ['Ez','Ex', 'Ey', 'Bz', 'Bx', 'By']):

    R1 = data['mesh']['R']
    Z1 = data['mesh']['Z']

    points1 = np.zeros([Z1.size, 2])
    points1[:, 0] = Z1.flatten()
    points1[:, 1] = R1.flatten()
    fields = dict()

    for field_tag in field_list:
        Z, R, F = getField(data, sol, field_tag)
        points0 = (Z[:, 0], R[0, :])
        fields[field_tag] = np.reshape(interpn(points0, F, points1, bounds_error=False, fill_value=None), np.shape(Z1))

    return fields

def plot(data, title = '', plot_type = 'default', scale = 'default', inPoints ='off', levels=100, clims=[0,100], saveFlag='off', **kwargs):

    if not('field' in kwargs):
        if data['geometry']['nz'] == data['mesh']['Z'].shape[0] - 1:  # New FEM meshes (not staggered):
            Z   = data['mesh']['Z']
            R   = data['mesh']['R']
            val = kwargs['solution'][kwargs['field_tag']]
        else:
            Z, R, val = getField(data, kwargs['solution'], kwargs['field_tag'])
    else:
        Z = kwargs['Z']
        R = kwargs['R']
        val = kwargs['field']

    if inPoints == 'on':
        if 'PointFlag' in data.keys():
            PointFlag = data['PointFlag']
        else:
            PointFlag = in_out(data, Z, R)
        valaux =val.flatten(); valaux[PointFlag == False] = float('nan')
        val = np.reshape(valaux, Z.shape)
        
    if plot_type == 'phase':
        val = np.angle(val) 
        title = title + ' Phase'

        plt.hsv()
        p = plt.contourf(Z, R, val, levels = levels)
        for c in p.collections:
            c.set_edgecolor("face")
        plt.title(title)
        plt.colorbar(p)

    else:
        if plot_type == 'magnitude':
            val = abs(val)
            title = title + ' Magnitude'  

        plt.viridis()
        if not(scale == 'log'):
            inter = np.nanmax(val) - np.nanmin(val)
            vmin  = np.float64(np.nanmin(val)+inter*clims[0]/100)
            vmax  = np.float64(np.nanmin(val)+inter*clims[1]/100)
            p = plt.contourf(Z, R, val, levels = levels, vmin=vmin, vmax=vmax)
        else:
            p = plt.contourf(Z, R, val, levels = levels, norm=LogNorm())
                
        for c in p.collections:
            c.set_edgecolor("face")
        plt.title(title)
        plt.colorbar(p)
    
    if (saveFlag=='on'): 
        if not(os.path.isdir(data['general']['simdir'] + '/figs')):
                os.mkdir(data['general']['simdir'] + '/figs')
        figure_dir = data['general']['simdir'] + '/figs/' + title + '_' + plot_type + '_' + scale + '.png'
        plt.savefig(figure_dir)
    # else:
    plt.draw()


def power_abs(data, sol, **kwargs):

    if not('Ez' in kwargs):
        fields = interpFields(data, sol, field_list = ['Ez'])
        Ez = fields['Ez']
    else:
        Ez = kwargs['Ez']
    if not('Ex' in kwargs):
        fields = interpFields(data, sol, field_list = ['Ex'])
        Ex = fields['Ex']
    else:
        Ex = kwargs['Ex']
    if not('Ey' in kwargs):
        fields = interpFields(data, sol, field_list = ['Ey'])
        Ey = fields['Ey']
    else:
        Ey = kwargs['Ey']

    # Compute power absorption
    eps0 = 8.854187817e-12
    omega = 2*np.pi*data['simulation']['freq']*1.0e6

    P = np.zeros(data['mesh']['Z'].shape)

    for iz in range(0, 2*data['geometry']['nz'] + 1):
        for ix in range(0, 2*data['geometry']['nx'] + 1):
            Evect = np.array([Ez[iz, ix], Ex[iz, ix], Ey[iz, ix]], dtype=complex)
            J = -1j*np.inner(data['kappa'][iz, ix, :, :] - np.eye(3), Evect)*eps0*omega
            P[iz, ix] = 0.5*np.real(np.inner(np.conj(J), Evect))
    return P

def in_out(data, Zq, Rq):

    z_offset  = 0.05

    dof_x = Zq.flatten()
    dof_y = Rq.flatten()

    filename = data['general']['simdir'] + '/CORE/out/SimState.hdf5'

    if os.path.isfile(filename):

        f = h5py.File(filename, 'r')
        vertices = np.array(f['picM']['points'])

        f.close()

        # Identify mesh points of W-module outside of hyphen simulation domain
        tol = 1e-5
        vertices[:, 0] += z_offset
        vertices[:, 1] -= tol

        p = path.Path(vertices)

        points = np.zeros([np.size(dof_x), 2])
        points[:, 0] = dof_x
        points[:, 1] = dof_y

        PointFlag = p.contains_points(points)
    else:
        PointFlag = np.full(dof_x.shape, True)

    return PointFlag

def export(data,sol):

    if data['geometry']['nz'] == data['mesh']['Z'].shape[0] - 1:
        exportFEM(data, sol, filename = data['general']['name'])
    else:
        exportFD(data, sol, filename = data['general']['name'])


def exportFD(data, sol, filename = '', field_list = ['omega_ce','omega_pe', 'nu', 'thetaB', 'Jz', 'Jx', 'Jy', 'Ez', 'Ex', 'Ey', 'Bz', 'Bx', 'By', 'P', 'CMA', 'time',
                                                   'pre_time','assembly_time','solve_time','PointFlag','nthreads'], exceptions = []):
    
    time          = data['general']['time']
    pre_time      = data['general']['pre_time']
    assembly_time = data['general']['assembly_time']
    solve_time    = data['general']['solve_time']

    omega_ce = data['plasma'][:, :, 0]
    omega_pe = data['plasma'][:, :, 2]
    nu       = data['plasma'][:, :, 4]
    thetaB   = data['plasma'][:, :, 6]

    Jz = data['current'][:, :, 0]
    Jx = data['current'][:, :, 1]
    Jy = data['current'][:, :, 2]

    if not('nthreads' in data['general']):
        nthreads = 0
    else:
        nthreads = data['general']['nthreads']

    for field_tag in exceptions:
        field_list.remove(field_tag)

    for field_tag in ['Ez', 'Ex', 'Ey', 'Bz', 'Bx', 'By']:
        if field_tag in field_list:
            Z, R, val = getField(data, sol, field_tag)
            locals()[field_tag] = val

    if 'CMA' in field_list:
        CMA = cma(data['mesh']['Z'], data['mesh']['R'], data['plasma'][:, :, 0], data['plasma'][:, :, 2], 0) #1/(1836*40)
    if 'PointFlag' in field_list:
        if 'PointFlag' in data.keys():
            PointFlag = data['PointFlag']
        else:
            PointFlag = in_out(data, data['mesh']['Z'], data['mesh']['R'])
    if 'P' in field_list:
        P = power_abs(data, sol)
            
    output = dict()
    output['Z'] = data['mesh']['Z']
    output['R'] = data['mesh']['R']
    for field_tag in field_list:
        output[field_tag] = locals()[field_tag]

    if not(filename == ''):
        while os.path.exists(data['general']['simdir'] + filename + '.h5'):
            filename += "_"
        hf = h5py.File(data['general']['simdir'] + filename + '.h5', 'w')
        for field_tag in output:
            hf.create_dataset(field_tag, data = output[field_tag])
        hf.close()

    return output


def exportFEM(data, fem_sol, filename = '', field_list = ['omega_ce','omega_pe', 'nu', 'thetaB', 'Jz', 'Jx', 'Jy', 'Ez', 'Ex', 'Ey', 'Qa', 'CMA', 'time',
                                                   'pre_time','assembly_time','solve_time','PointFlag','nthreads'], exceptions = []):
    
    time          = data['general']['time']
    pre_time      = data['general']['pre_time']
    assembly_time = data['general']['assembly_time']
    solve_time    = data['general']['solve_time']

    omega_ce = data['plasma'][:, :, 0]
    omega_pe = data['plasma'][:, :, 2]
    nu       = data['plasma'][:, :, 4]
    thetaB   = data['plasma'][:, :, 6]

    Ez = fem_sol['Ez']
    Ex = fem_sol['Ex']
    Ey = fem_sol['Ey']
    Qa = fem_sol['Qa']

    Jz = data['current'][:, :, 0]
    Jx = data['current'][:, :, 1]
    Jy = data['current'][:, :, 2]

    if not('nthreads' in data['general']):
        nthreads = 0
    else:
        nthreads = data['general']['nthreads']

    for field_tag in exceptions:
        field_list.remove(field_tag)

    if 'CMA' in field_list:
        CMA = cma(data['mesh']['Z'], data['mesh']['R'], data['plasma'][:, :, 0], data['plasma'][:, :, 2], 0) #1/(1836*40)
    if 'PointFlag' in field_list:
        if 'PointFlag' in data.keys():
            PointFlag = data['PointFlag']
        else:
            PointFlag = in_out(data, data['mesh']['Z'], data['mesh']['R'])
            
    output = dict()
    output['Z'] = data['mesh']['Z']
    output['R'] = data['mesh']['R']
    for field_tag in field_list:
        output[field_tag] = locals()[field_tag]

    if not(filename == ''):
        while os.path.exists(data['general']['simdir'] + filename + '.h5'):
            filename += "_"
        hf = h5py.File(data['general']['simdir'] + filename + '.h5', 'w')
        for field_tag in output:
            hf.create_dataset(field_tag, data = output[field_tag])
        hf.close()

    return output

def cma(Zm, Rm, omega_ce, omega_pe, lamb):

    omega_pi = omega_pe*np.sqrt(lamb)
    omega_ci = omega_ce*lamb

    R = 1 - omega_pe**2/(1 - omega_ce) - omega_pi**2/(1 + omega_ci)
    L = 1 - omega_pe**2/(1 + omega_ce) - omega_pi**2/(1 - omega_ci)
    P = 1 - omega_pe**2 - omega_pi**2

    # S = 1 - omega_pe**2/(1-omega_ce**2)
    # D = - omega_pe**2*omega_ce/(1 - omega_ce**2)

    S = 1/2*(R + L)
    D = 1/2*(R - L)

    log = np.zeros([Zm.shape[0], Zm.shape[1], 8])
    log[:, :, 6] = np.ones([Zm.shape[0], Zm.shape[1]])
    
    # Region 1
    log[:, :, 0] = R > 0
    log[:, :, 1] = L > 0
    log[:, :, 2] = P > 0
    log[:, :, 3] = S >= 0
    log[:, :, 4] = D <= 0
    log[:, :, 5] = P/S > 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond1 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))
    # cond1 = D < 0

    # Region 2
    log[:, :, 0] = R < 0
    log[:, :, 1] = L > 0
    log[:, :, 2] = P > 0
    log[:, :, 3] = S > 0
    log[:, :, 4] = D < 0
    log[:, :, 5] = P/S > 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond2 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    # Region 3
    log[:, :, 0] = R < 0
    log[:, :, 1] = L > 0
    log[:, :, 2] = P > 0
    log[:, :, 3] = S < 0
    log[:, :, 4] = D < 0
    log[:, :, 5] = P/S < 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond3 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    # Region 4
    log[:, :, 0] = R < 0
    log[:, :, 1] = L > 0
    log[:, :, 2] = P < 0
    log[:, :, 3] = S < 0
    log[:, :, 4] = D < 0
    log[:, :, 5] = P/S > 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond4 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    # Region 5
    log[:, :, 0] = R < 0
    log[:, :, 1] = L < 0
    log[:, :, 2] = P < 0
    log[:, :, 3] = S <= 0
    log[:, :, 4] = D <= 0
    log[:, :, 5] = P/S > 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond5 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    # Region 6
    log[:, :, 0] = R > 0
    log[:, :, 1] = L > 0
    log[:, :, 2] = P > 0
    log[:, :, 3] = S > 0
    log[:, :, 4] = D > 0
    log[:, :, 5] = P/S > 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond6 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    # Region 7
    log[:, :, 0] = R > 0
    log[:, :, 1] = L > 0
    log[:, :, 2] = P < 0
    log[:, :, 3] = S > 0
    log[:, :, 4] = D > 0
    log[:, :, 5] = P/S < 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond7 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    # Region 8
    log[:, :, 0] = R > 0
    log[:, :, 1] = L < 0
    log[:, :, 2] = P < 0
    log[:, :, 3] = S > 0
    log[:, :, 4] = D > 0
    log[:, :, 5] = P/S < 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond8 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    F = np.zeros(Zm.shape)

    for iz in range(0, Zm[:, 0].size):
        for ix in range(0, Zm[0, :].size):
            if cond1[iz, ix] == 1:
                F[iz, ix] = 1
            elif cond2[iz, ix] == 1:
                F[iz, ix] = 2
            elif cond3[iz, ix] == 1:
                F[iz, ix] = 3
            elif cond4[iz, ix] == 1:
                F[iz, ix] = 4
            elif cond5[iz, ix] == 1:
                F[iz, ix] = 5
            elif cond6[iz, ix] == 1:
                F[iz, ix] = 6
            elif cond7[iz, ix] == 1:
                F[iz, ix] = 7
            elif cond8[iz, ix] == 1:
                F[iz, ix] = 8
            else:
                F[iz, ix] = np.nan

    # scp.gaussian_filter(F , 10)
    return F