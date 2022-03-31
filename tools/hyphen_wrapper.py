"""
Created on Fri Aug  9 20:07:51 2019
@author: jiewei
"""

import numpy as np 
import h5py
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import NearestNDInterpolator
from matplotlib import path
import scipy.ndimage.filters as scp
import math


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
    # Initialize output variable
    Nz = 2*data['geometry']['nz'] + 1
    Nx = 2*data['geometry']['nx'] + 1
    plasma = np.zeros([Nz, Nx, 9], dtype = complex)
    
    # Prepare W-module mesh

    R = data['mesh']['R']
    Z = data['mesh']['Z']
    
    z_offset = data['geometry']['z_offset']
    dof_x = Z.flatten()
    dof_y = R.flatten()

    # Load data from hyphen
    r = open(sim_path + 'CORE/inp/' + 'sim_params.inp','r')
    lines = r.readlines()
    r.close() 
    
    nlines = len(lines) 
    for i in range(0,nlines):
        if "B_fact" in lines[i]:
            B_fact = float(lines[i][lines[i].find('=')+1:-1])  
    
    
    
    filename = sim_path + 'CORE/out/' + 'SimState.hdf5'
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

    fNu  = LinearNDInterpolator(coord_fld, Nu)
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

    ## Obtain magnetic field
    BZwave, BRwave, normB = get_B(sim_path, data['simulation']['B_file_name'], z_offset, dof_x, dof_y)
    BZwave = B_fact*BZwave
    BRwave = B_fact*BRwave
    normB  = B_fact*normB
    
    BZwave[np.isnan(BZwave)] = 0
    BRwave[np.isnan(BRwave)] = 0
    normB[np.isnan(normB)]   = 0

    # Obtain non-dimensional properties
    e_charge = 1.60217662e-19
    e_mass   = 9.10938356e-31
    eps0     = 8.854187817e-12
    
    omega    = 2*np.pi*data['simulation']['freq']*1.0e6
    omega_pe = np.reshape(np.sqrt(Nwave*e_charge**2/e_mass/eps0)/omega,(Nz, Nx))
    omega_ce = np.reshape((e_charge*normB/e_mass)/omega,(Nz, Nx))
    nu       = np.reshape(NUwave/omega,(Nz, Nx))
    
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
    
    filename = sim_path + 'CORE/out/' + 'SimState.hdf5'
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
    

    # Open the magnetic field file
    r = open(sim_path + 'SET/inp/' + B_file_name, "r")
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
       
    fBZ = LinearNDInterpolator(points, Bz_points, fill_value = 0)    
    fBR = LinearNDInterpolator(points, Br_points, fill_value = 0)    
    
    Bz = fBZ(zs, rs) 
    Br = fBR(zs, rs)
  
    B = np.sqrt(Bz**2.0 + Br**2.0)
    
    
    return Bz, Br, B