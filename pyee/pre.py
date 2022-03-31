import sys
import yaml
import numpy as np
import pickle

def prepareData(path, mode, user_data):
    with open(path + "simrc.yml", 'r') as stream:
        data = yaml.safe_load(stream)

    # Modify selected data
    for entry in data:
        if entry in user_data and type(data[entry] == 'dict'):
            for entry1 in data[entry]:
                if entry1 in user_data[entry]: 
                    data[entry][entry1] = user_data[entry][entry1]
        elif entry in user_data:
            data[entry] = user_data[entry]

    if mode == 'solve_pre':  # FEM pre can not be load by FD solver (the other way around works)
        print("Loading preprocessor file, configuration will be overwritten")
        general_temp = data['general'].copy()                                           # I don't want this overwritten
        file = open(data['general']['simdir'] + "/preproc.pkl","rb")
        data = pickle.load(file)
        data['general'] = general_temp
        file.close()
    else:
        # Create mesh
        data['mesh'] = dict()

        if data['general']['fem'] == True:
            Nz = data['geometry']['nz'] + 1
            Nx = data['geometry']['nx'] + 1
        else:
            Nz = 2*data['geometry']['nz'] + 1
            Nx = 2*data['geometry']['nx'] + 1

        R1, Z1 = np.meshgrid(np.linspace(0, data['geometry']['Lx'], Nx), np.linspace(0, data['geometry']['Lz'],Nz))

        data['mesh']['R'] = R1
        data['mesh']['Z'] = Z1

        data['plasma'] = np.zeros([Nz, Nx, 8], dtype = complex)

        if data['conditions']['kappa'][0] == 'file':
            filename = data['conditions']['kappa'][1]
            sys.path.append(filename + '/')
            import config_sim
            data['plasma']  = config_sim.plasma(data)
            sys.path.remove(filename + '/')
            del sys.modules['config_sim']
        data['kappa']   = kappa(data)

        if data['conditions']['current'][0] == 'file':
            filename = data['conditions']['current'][1]
            sys.path.append(filename + '/')
            import config_sim
            data = config_sim.current(data)
            sys.path.remove(filename + '/')
            del sys.modules['config_sim']
        elif data['conditions']['current'][0] == 'function':
            data['current'] = current(data)
        else:
            raise Exception('Not implemented')

        data['BCs'] = boundaries(data)

        file = open(data['general']['simdir'] + "/preproc.pkl","wb")
        pickle.dump(data,file)
        file.close()

    return data

def kappa(data):
    params = data['conditions']['kappa']
    
    dx = data['geometry']['Lx'] / data['geometry']['nx']
    dz = data['geometry']['Lz'] / data['geometry']['nz']

    if data['geometry']['nz'] == data['mesh']['Z'].shape[0] - 1:  
        rangez = range(0, data['geometry']['nz'] + 1)
        rangex = range(0, data['geometry']['nx'] + 1)
        kappa  = np.zeros([data['geometry']['nz'] + 1, data['geometry']['nx'] + 1, 3, 3], dtype = complex)
    else:                                                 # Legacy FD staggered
        rangez = range(0, 2*data['geometry']['nz'] + 1)
        rangex = range(0, 2*data['geometry']['nx'] + 1)
        kappa  = np.zeros([2*data['geometry']['nz'] + 1, 2*data['geometry']['nx'] + 1, 3, 3], dtype = complex)

    for iz in rangez:
        for ix in rangex:

            z = data['mesh']['Z'][iz,ix]
            x = data['mesh']['R'][iz,ix]

            if params[0] == 'function':
                omega_ce = abs(eval(params[1]))  # Force unsigned
                omega_ci = abs(eval(params[2]))
                omega_pe = eval(params[3])
                omega_pi = eval(params[4])
                nu_e     = eval(params[5])
                nu_i     = eval(params[6])
                thetaB   = eval(params[7])
                phiB     = eval(params[8])
            elif params[0] == 'file':  
                omega_ce = data['plasma'][iz, ix, 0]
                omega_ci = data['plasma'][iz, ix, 1]
                omega_pe = data['plasma'][iz, ix, 2]
                omega_pi = data['plasma'][iz, ix, 3]
                nu_e     = data['plasma'][iz, ix, 4]
                nu_i     = data['plasma'][iz, ix, 5]
                thetaB   = data['plasma'][iz, ix, 6]
                phiB     = data['plasma'][iz, ix, 7]
            else:
                raise Exception('Not implemented')

            # Stix components
            R = 1 - omega_pe**2 / (1 + 1j * nu_e - omega_ce) - omega_pi**2 / (1 + 1j * nu_i + omega_ci)
            L = 1 - omega_pe**2 / (1 + 1j * nu_e + omega_ce) - omega_pi**2 / (1 + 1j * nu_i - omega_ci)
            P = 1 - omega_pe**2 / (1 + 1j * nu_e) - omega_pi**2 / (1 + 1j * nu_i)
            S = (R + L) / 2
            D = (R - L) / 2

            # Kappa tensor in Stix reference frame
            k = np.array([[P, 0, 0], [0, S, -1j * D], [0, 1j * D, S]], dtype = complex)

            # Rotation matrices
            Rtheta = np.array([[np.cos(thetaB), np.sin(thetaB), 0], [-np.sin(thetaB), np.cos(thetaB), 0], [0, 0, 1]], dtype = complex)
            Rphi   = np.array([[np.cos(phiB), 0, np.sin(phiB)], [0, 1, 0], [-np.sin(phiB), 0, np.cos(phiB)]], dtype = complex)

            # Rotated tensor
            left = np.matmul(np.transpose(Rtheta),np.transpose(Rphi))
            right = np.matmul(Rphi, Rtheta)
            kappa[iz, ix, :, :] = np.matmul(left, np.matmul(k, right))

    return kappa

def current(data):

    params = data['conditions']['current']

    if data['geometry']['nz'] == data['mesh']['Z'].shape[0] - 1:  
        rangez = range(0, data['geometry']['nz'] + 1)
        rangex = range(0, data['geometry']['nx'] + 1)
        ja     = np.zeros([data['geometry']['nz'] + 1, data['geometry']['nx'] + 1, 3], dtype = complex)
    else:                                                 # Legacy FD staggered
        rangez = range(0, 2*data['geometry']['nz'] + 1)
        rangex = range(0, 2*data['geometry']['nx'] + 1)
        ja     = np.zeros([2*data['geometry']['nz'] + 1, 2*data['geometry']['nx'] + 1, 3], dtype = complex)


    if params[0] == 'function':
        for iz in rangez:
            for ix in rangex:

                z = data['mesh']['Z'][iz,ix]
                x = data['mesh']['R'][iz,ix]

                ja[iz, ix, 0] = eval(params[1])
                ja[iz, ix, 1] = eval(params[2])
                ja[iz, ix, 2] = eval(params[3])
    else:
        raise Exception('Not implemented')
    return ja

def boundaries(data):

    Lz = data['geometry']['Lz']
    Lx = data['geometry']['Lx']
    nz = data['geometry']['nz']
    nx = data['geometry']['nx']
    dz = Lz / nz
    dx = Lx / nx

    EqEz = np.ones((nz + 1,nx + 1),np.int8)
    EqEx = np.ones((nz + 1,nx + 1),np.int8)
    EqEy = np.ones((nz + 1,nx + 1),np.int8)

    if not('boundaries' in data):
        num_rect = 0
    else:
        num_rect = len(data['boundaries']['zo'])

    for j in range(0, num_rect):

        zo       = data['boundaries']['zo'][j]
        ro       = data['boundaries']['xo'][j]
        L        = data['boundaries']['w'][j]
        H        = data['boundaries']['h'][j]

        iz0 = round(zo/dz)
        ir0 = round(ro/dx)
        nL  = round(L/dz)
        nH  = round(H/dx) 

        for i in range(nL): #HORIZONTAL
            if not(ir0 == 0): # Don't allow horizontal walls at the axis
                EqEz[iz0 + i, ir0] = 0
                EqEy[iz0 + i, ir0] = 0
                EqEy[iz0 + i + 1, ir0] = 0

            EqEz[iz0 + i, ir0 + nH] = 0
            EqEy[iz0 + i, ir0 + nH] = 0
            EqEy[iz0 + i + 1, ir0 + nH] = 0

        for i in range(nH):
            EqEx[iz0, ir0 + i] = 0
            EqEy[iz0, ir0 + i] = 0
            EqEy[iz0, ir0 + i + 1] = 0

            EqEx[iz0 + nL, ir0 + i] = 0
            EqEy[iz0 + nL, ir0 + i] = 0
            EqEy[iz0 + nL, ir0 + i + 1] = 0


    BCs = dict()

    BCs['EqEz'] = EqEz
    BCs['EqEx'] = EqEx
    BCs['EqEy'] = EqEy

    return BCs