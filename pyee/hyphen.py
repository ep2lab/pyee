
import numpy as np 
import h5py
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import interpn
from matplotlib import path
import scipy.ndimage.filters as scp
import math
import yaml
import pyee
import os
import pickle

def pre(path):

    """

    ###############################################################################
    Description:    Read the input data from user and hyphen    
    ###############################################################################
    Inputs:         1) path: path of the simulation folder
    ###############################################################################
    Outputs:        1) data: input data

    
    """
    
    # Read user data
    with open(path + 'simrc.yml', 'r') as stream:
        data = yaml.safe_load(stream)

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

    # Interpolate plasma transport properties 
    data['plasma']  = interpH2W(data,path)
    
    # Compute dielectric tensor
    data['conditions'] = dict()
    data['conditions']['kappa'] = dict()
    data['conditions']['kappa'][0] = 'file'
    
    data['kappa']   = pyee.pre.kappa(data)

    # Define the antenna current density
    import tools

    antenna_data = dict()
    antenna_data['za'] = 0.075 + 0.2
    antenna_data['La'] = 0.075
    antenna_data['ra'] = 0.0175
    antenna_data['sa'] = 0.003
    antenna_data['h']  = 0.5      

    ja  = tools.antenna.helix(data, antenna_data, data['simulation']['modes'][0])
    data['current'] = ja

    data['BCs'] = pyee.pre.boundaries(data)

    file = open(data['general']['simdir'] + "/preproc.pkl","wb")
    pickle.dump(data,file)
    file.close()

    return data


def post(path,data,v):

    """

    ###############################################################################
    Description:    Save the pyee solution and update power depostion in hyphen   
    ###############################################################################
    Inputs:         1) path: path of the simulation folder
    Inputs:         2) data: input data
    Inputs:         3) v: solution from pyee
    
    """
    if data['geometry']['nz'] == data['mesh']['Z'].shape[0] - 1: # Check FEM or FD mesh
        Ez = v['Ez']
        Ex = v['Ex']
        Ey = v['Ey']
        Qa = np.abs(v['Qa'])
    else:
        Ex,Ey,Ez,Qa = extract_sol(data,v)
    
    # Obtain power absorption in E-module 
    P_EM,P_abs = interpW2H(path,data,Qa)
    
    # Scale the power 
    p_fact = data['simulation']['P_inp']/P_abs
    
    Qa = p_fact*Qa
    Qa_EM = p_fact*P_EM
    
    Ex = math.sqrt(p_fact)*Ex
    Ey = math.sqrt(p_fact)*Ey
    Ez = math.sqrt(p_fact)*Ez
    
    # Save solution
    R = data['mesh']['R']
    Z = data['mesh']['Z']
            
    h5W = h5py.File(path+'PostDataW.hdf5','w')
    
    h5W.create_dataset("/mesh/R",R.shape, h5py.h5t.NATIVE_DOUBLE)
    h5W.create_dataset("/mesh/Z",Z.shape, h5py.h5t.NATIVE_DOUBLE)
    h5W.create_dataset("/Ex",Ex.shape, dtype ='complex128')
    h5W.create_dataset("/Ey",Ey.shape, dtype ='complex128')
    h5W.create_dataset("/Ez",Ez.shape, dtype ='complex128')
    h5W.create_dataset("/Qa",Qa.shape, h5py.h5t.NATIVE_DOUBLE)
    
    dataset = h5W["/mesh/R"]
    dataset[...] = R
    dataset = h5W["/mesh/Z"]
    dataset[...] = Z
    dataset = h5W["/Ex"]
    dataset[...] = Ex
    dataset = h5W["/Ey"]
    dataset[...] = Ey
    dataset = h5W["/Ez"]
    dataset[...] = Ez
    dataset = h5W["/Qa"]
    dataset[...] = Qa
    
    
    h5W.close()   
    
    
    h5H = h5py.File(path+'SimState.hdf5','r+')
    
    dataset = h5H["/ssD_eFld_e_inst/Qa"]
    dataset[...] = Qa_EM
    
    
    h5H.close()  


def interpH2W(data,sim_path):
    
    """
    
    ###############################################################################
    Description:    Interpolate necessary plasma data from hyphen to pyee
    ###############################################################################
    Inputs:         1) data: input data
    Inputs:         2) sim_path: path of the simulation folder
    ###############################################################################
    Outputs:        1) plasma: plasma transport properties from hyphen
        
    """    
    # Prepare W-module mesh
    R = data['mesh']['R']
    Z = data['mesh']['Z']

    # Initialize output variable
    Nz = Z.shape[0]
    Nx = Z.shape[1]
    plasma = np.zeros([Nz, Nx, 9], dtype = complex)
    
    z_offset = data['geometry']['z_offset']
    dof_x = Z.flatten()
    dof_y = R.flatten()

    # Load data from hyphen
    r = open(sim_path + 'sim_params.inp','r')
    lines = r.readlines()
    r.close() 
    
    nlines = len(lines) 
    for i in range(0,nlines):
        if "B_fact" in lines[i]:
            B_fact = float(lines[i][lines[i].find('=')+1:-1])  
    
    filename = sim_path + 'SimState.hdf5'

    f = h5py.File(filename, 'r')
    
    ZPic = np.array(f['picM']['zs']).flatten() + z_offset
    RPic = np.array(f['picM']['rs']).flatten()
    N = np.array(f['ssD_picM_acc']['n']).flatten()

    Zfld = np.concatenate((np.array(f['eFldM']['element_geom'])[0, :],
           np.array(f['eFldM']['face_geom'])[0, :]), 0) + z_offset     
    Rfld = np.concatenate((np.array(f['eFldM']['element_geom'])[1, :],
           np.array(f['eFldM']['face_geom'])[1, :]), 0)
    Nu   = np.concatenate((np.array(f['ssD_eFld_e_acc']['freq_e_tot']),
                          np.array(f['ssD_eFld_f_acc']['freq_e_tot'])), 1).flatten()
    
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

    # Define interpolators
    coord = np.zeros([np.size(ZPic),2])
    coord[:, 0] = ZPic
    coord[:, 1] = RPic

    fN  = LinearNDInterpolator(coord, N, fill_value = 0)


    coord_fld = np.zeros([np.size(Zfld),2])
    coord_fld[:, 0] = Zfld
    coord_fld[:, 1] = Rfld

    fNu = LinearNDInterpolator(coord_fld, Nu)
    fNuN = NearestNDInterpolator(coord_fld, Nu)

    Nwave  = fN(dof_x, dof_y)
    NUwave = fNu(dof_x, dof_y)
    
    
    PointFlagRS = np.reshape(PointFlag,(Nz, Nx))
    for iz in range(0,Nz):
        for ir in range(0,Nx):
            if PointFlagRS[iz,ir] == True:
                if math.isnan(NUwave[iz*Nx+ir])==True:
                    NUwave[iz*Nx+ir] = fNuN(dof_x[iz*Nx+ir], dof_y[iz*Nx+ir])
    
    Nwave[PointFlag == False]  = 0
    NUwave[PointFlag == False] = 0
        

    BZwave, BRwave, normB = get_B(sim_path, "B.dat", z_offset, dof_x, dof_y)
    BZwave = B_fact*BZwave
    BRwave = B_fact*BRwave
    normB  = B_fact*normB

    # Obtain non-dimensional properties
    e_charge = 1.60217662e-19
    e_mass   = 9.10938356e-31
    eps0     = 8.854187817e-12
    
    omega     = 2*np.pi*data['simulation']['freq']*1.0e6
    omega_pe  = np.reshape(np.sqrt(Nwave*e_charge**2/e_mass/eps0)/omega,(Nz, Nx))
    omega_ce  = np.reshape((e_charge*normB/e_mass)/omega,(Nz, Nx))
    nu        = np.reshape(NUwave/omega,(Nz, Nx))
    PointFlag = np.reshape(PointFlag, Z.shape)

    # Filter  for a smooth transition
    min_omg_pe = np.sqrt(0e14*e_charge**2/eps0/e_mass)/omega
    f = 0
    omega_pe[omega_pe < min_omg_pe]  = min_omg_pe
    omega_pe = scp.gaussian_filter(omega_pe, f*np.array(Z.shape)/100)
    # omega_ce = scp.gaussian_filter(omega_ce, f*np.array(Zw.shape)/100)

    # Artificial collision band
    # nu[np.logical_and(omega_pe>0.75, omega_pe<1.5)] = 100
    # nu = scp.gaussian_filter(nu, (1*Nz/100, 1*Nx/100))

    plasma[:, :, 0] = omega_ce
    plasma[:, :, 1] = 0
    plasma[:, :, 2] = omega_pe
    plasma[:, :, 3] = 0
    plasma[:, :, 4] = nu
    plasma[:, :, 5] = 0
    plasma[:, :, 6] = np.reshape(np.arctan2(BRwave, BZwave),(Nz, Nx))
    plasma[:, :, 7] = 0
    plasma[:, :, 8] = PointFlagRS

    return plasma


def interpW2H(sim_path,data,P_WM):
    
    """
    
    ###############################################################################
    Description:    Interpolate the power absorption from pyee to hyphen
    ###############################################################################
    Inputs:         1) sim_path: path of the simulation folder
    Inputs:         2) data: input data    
    Inputs:         3) P_WM: power absorption in pyee mesh
    ###############################################################################
    Outputs:        1) P_EM: power absorption in E-module mesh
        
    """    
    
    # Load E-module mesh
    z_offset = data['geometry']['z_offset']
    
    filename = sim_path + 'SimState.hdf5'
    f = h5py.File(filename, 'r')
    
    Zfld_e = np.array(f['eFldM']['element_geom'])[0, :] + z_offset     
    Rfld_e = np.array(f['eFldM']['element_geom'])[1, :]
    
    Vol_e = np.array(f['eFldM']['element_geom'])[4, :]

    
    f.close()
    
    
    # Initialize output variable
    P_EM = np.zeros([np.size(Zfld_e),1])
    
    
    # Prepare W-module mesh

    R = data['mesh']['R']
    Z = data['mesh']['Z']
    
   
    dof_x = Z.flatten()
    dof_y = R.flatten()
    
    # Define interpolator
    
    coord = np.zeros([np.size(dof_x),2])
    coord[:, 0] = dof_x
    coord[:, 1] = dof_y
    
    p_wm = P_WM.flatten()

    fP = LinearNDInterpolator(coord, p_wm, fill_value = 0)
    
    P_EM  = fP(Zfld_e, Rfld_e)
    
    P_abs = 0.0

 
    for i in range(0,np.size(Zfld_e)):
        P_abs = P_abs+P_EM[i]*Vol_e[i]
   
             

    return P_EM,P_abs

def get_B(sim_path, B_file_name, offset_zB, zs, rs):
    
    """

    ###############################################################################
    Description:    Read and interpolate the magnetic field    
    ###############################################################################
    Inputs:         1) sim_path: path of the simulation
                    2) B_file_name: name of the input B file with extension
                    3) offset_zB: offset on z coordinate of HDF5 B file
                    4) z,r: meshgrid matrices
    
    ###############################################################################
    Outputs:        1) Bz, Br: z and r components of the B field at the PIC mesh nodes
                    2) B: B field magnitude at the PIC mesh nodes   

    
    """

    import numpy as np
    from scipy.interpolate import LinearNDInterpolator
    

    if os.path.isfile(sim_path + "/Bpre.pkl"):
        file = open(sim_path + "/Bpre.pkl","rb")
        Bdict = pickle.load(file)
        file.close()
    else:
        print('Creating interpolators for the magnetic field. The first call to the W-module can be slow')
        # Open the magnetic field file
        r = open(sim_path + B_file_name, "r")
        r.readline() 
        firstline = r.readline()
        pos1 = firstline.find('[')
        pos2 = firstline.find(']')
        pos3 = pos2 + 1 + firstline[pos2+1:].find('[')
        pos4 = pos2 + 1 + firstline[pos2+1:].find(']')
        pos5 = pos4 + 1 + firstline[pos4+1:].find('[')
        pos6 = pos4 + 1 + firstline[pos4+1:].find(']')
        points_min = firstline[pos1+1:pos2].split(" ")
        rmin = float(points_min[0])
        zmin = float(points_min[1])
        points_max = firstline[pos3+1:pos4].split(" ")
        rmax = float(points_max[0])
        zmax = float(points_max[1])
        deltas = firstline[pos5+1:pos6].split(" ")
        delta_r = float(deltas[0])
        delta_z = float(deltas[1])
        lines = r.readlines()
        r.close()
            
        npoints_r = int((rmax - rmin) / delta_r)+1
        npoints_z = int((zmax - zmin) / delta_z)+1
            
        points    = np.zeros((int(npoints_r*npoints_z),2),dtype='float')
        Br_points = np.zeros((int(npoints_r*npoints_z)),dtype='float')
        Bz_points = np.zeros((int(npoints_r*npoints_z)),dtype='float')
        
        # Read and interpolate   
        for i in range(0,int(npoints_r*npoints_z)-1):
            data_list = lines[i][0:lines[i].find('\t')].split(" ")
            points[i,0] = float(data_list[1])+offset_zB
            points[i,1] = float(data_list[0])
            Bz_points[i] = float(data_list[5])    
            Br_points[i] = float(data_list[4])
        
        Bdict = dict()
        Bdict['Z'] = LinearNDInterpolator(points, Bz_points, fill_value = 0)    
        Bdict['R'] = LinearNDInterpolator(points, Br_points, fill_value = 0)    

        file = open(sim_path + "/Bpre.pkl","wb")      # Save interpolators for use in another time step
        pickle.dump(Bdict,file)
        file.close()
    
    Bz = Bdict['Z'](zs, rs) 
    Br = Bdict['R'](zs, rs)
  
    B = np.sqrt(Bz**2.0 + Br**2.0)
    
    return Bz, Br, B


def extract_sol(data, sol):
    
    """

    ###############################################################################
    Description:    Extract the pyee solution  (Only FD)
    ###############################################################################
    Inputs:         1) data: input data
    Inputs:         2) sol: solution from pyee
    ###############################################################################
    Outputs:        1) Ex
    Outputs:        2) Ey
    Outputs:        3) Ez
    Outputs:        4) P: power density deposited

    
    """
    # Extract the electric field and interpolate
    R = data['mesh']['R']
    Z = data['mesh']['Z']

    points = np.zeros([Z.size, 2])
    points[:, 0] = Z.flatten()
    points[:, 1] = R.flatten()

    Z3, R3, Ez = pyee.post.getField(data, sol, "Ez")
    points3 = (Z3[:, 0], R3[0, :])  
    Ez = np.reshape(interpn(points3, Ez, points, bounds_error=False), np.shape(Z))
 

    Z1, R1, Ex = pyee.post.getField(data, sol, "Ex")
    points1 = (Z1[:, 0], R1[0, :])
    Ex = np.reshape(interpn(points1, Ex, points, bounds_error=False), np.shape(Z))
  

    Z2, R2, Ey = pyee.post.getField(data, sol, "Ey")
    points2 = (Z2[:, 0], R2[0, :])
    Ey = np.reshape(interpn(points2, Ey, points, bounds_error=False), np.shape(Z))
 

    # Compute power absorption
    #eps0 = 8.854187817e-12
    #omega = 2*np.pi*data['simulation']['freq']

    P = np.zeros(Z.shape)

    for iz in range(0, 2*data['geometry']['nz'] + 1):
        for ix in range(0, 2*data['geometry']['nx'] + 1):
            Evect = np.array([Ez[iz, ix], Ex[iz, ix], Ey[iz, ix]], dtype=complex)
            J = -1j*np.inner(data['kappa'][iz, ix, :, :] - np.eye(3), Evect)#*omega*eps0
            P[iz, ix] = 0.5*np.real(np.inner(np.conj(J), Evect))
    return Ex,Ey,Ez,P