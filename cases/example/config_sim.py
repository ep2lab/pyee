import numpy as np 

# This optional file is used for complex configurations such as interpolating from external files,
# using complex antennas etc..
# Its usage can be activated in simrc.yml
# You can import any functions needed, filters etc..

def current(data):

    # Example for a half turn helical antenna
    import tools.antenna as antenna

    antenna_data = dict()
    antenna_data['za'] = 0.3       # Axial position of the center of the antenna
    antenna_data['La'] = 0.1       # Antenna axial length
    antenna_data['ra'] = 0.5       # Radial position of the antenna
    antenna_data['sa'] = 0.05      # Antenna current radial width
    antenna_data['h']  = 0.5       # Helix number

    ja  = antenna.helix(data, antenna_data, data['simulation']['modes'][0])

    data['current'] = ja

    return data

def plasma(data):

    R = data['mesh']['R']
    Z = data['mesh']['Z']

    # Initialize output variable
    Nz = Z.shape[0]
    Nx = Z.shape[1]
    plasma = np.zeros([Nz, Nx, 9], dtype = complex)

    # This dummy example is just a perfect vacuum
    plasma[:, :, 0] = 0         # Electron Cyclotron Frequency
    plasma[:, :, 1] = 0         # Ion Cyclotron Frequency
    plasma[:, :, 2] = 0         # Electron Plasma Frequency
    plasma[:, :, 3] = 0         # Ion Plasma Frequency
    plasma[:, :, 4] = 0         # Effective Electron Collision Frequency
    plasma[:, :, 5] = 0         # Effective Ion Collision Frequency
    plasma[:, :, 6] = 0         # Background Magnetic Field Angle in the z-x plane
    plasma[:, :, 7] = 0         # Background Magnetic Field Phi Angle (spherical, rarely used)

    return plasma