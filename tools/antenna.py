import numpy as np 

def helix(data, antenna_data, mode):

    if data['geometry']['nz'] == data['mesh']['Z'].shape[0] - 1:
        rangez = range(1, data['geometry']['nz'])
        rangex = range(1, data['geometry']['nx'])
        ja     = np.zeros([data['geometry']['nz'] + 1, data['geometry']['nx'] + 1, 3], dtype = complex)
    else:                                                 # Legacy FD staggered
        rangez = range(1, 2*data['geometry']['nz'])
        rangex = range(1, 2*data['geometry']['nx'])
        ja  = np.zeros([2*data['geometry']['nz'] + 1, 2*data['geometry']['nx'] + 1, 3], dtype = complex)
    
    m = mode
    
    # Antenna parameters
    za = antenna_data['za']
    La = antenna_data['La']
    ra = antenna_data['ra']
    sa = antenna_data['sa']
    h  = antenna_data['h']
              
    rsamp = np.linspace(0, data['geometry']['Lx'], num = 2 * data['geometry']['nx'] + 1)
    fac = np.trapz(np.exp(-(rsamp - ra)**2/sa**2), x = rsamp)

    for iz in rangez:
        for ix in rangex:

            z = data['mesh']['Z'][iz, ix]
            x = data['mesh']['R'][iz, ix]

            if z >= (za - La/2) and z <= (za + La/2):
                ja[iz, ix, 0] = (np.exp(-(x - ra)**2/sa**2)/fac)*(1-np.exp(-1j*m*np.pi))/x*((2*np.pi*h*x/La)**2+1)**(-0.5)*np.exp(-2*1j*m*np.pi*h*(z-za+La/2)/La)/(2*np.pi)

    for iz in rangez:
        for ix in rangex:
            x    = data['mesh']['R'][iz, ix]
            dz_2 = data['mesh']['Z'][iz + 1, ix] - data['mesh']['Z'][iz - 1, ix]
            ja[iz, ix, 2] = 1j*x*(ja[iz + 1, ix, 0] - ja[iz - 1, ix, 0])/(dz_2)/m

    # if data['geometry']['nz'] == data['mesh']['Z'].shape[0] - 1:  # New FEM meshes (not staggered)
    #     ja = ja[0::2,0::2,:]

    return ja