'''
Creation Date: Wednesday March 4th 2020
Author: Pedro Jimenez, pejimene@ing.uc3m.es
##################################################
File: test.py
Project: pyee
##################################################
Last Modified: Monday, April 20th 2020, 6:13:10 pm
Modified By: Pedro Jimenez
##################################################
Code Description:
##################################################
Copyright (c) 2020 Plasma & Space Propulsion Team (EP2). Universidad Carlos III de Madrid.
'''
from scipy.interpolate import interpn
import numpy as np
import run_pyee
import pyee

data, sol = run_pyee.sim('helicon')

R1 = data['mesh']['R']
Z1 = data['mesh']['Z']

points1 = np.zeros([Z1.size, 2])
points1[:, 0] = Z1.flatten()
points1[:, 1] = R1.flatten()

Z, R, Ez = pyee.post.getField(data, sol, "Ez")
points0 = (Z[:, 0], R[0, :])
Ez = np.reshape(interpn(points0, Ez, points1, bounds_error=False, fill_value=0), np.shape(Z1))

Z, R, Ex = pyee.post.getField(data, sol, "Ex")
points0 = (Z[:, 0], R[0, :])
Ex = np.reshape(interpn(points0, Ex, points1, bounds_error=False, fill_value=0), np.shape(Z1))

Z, R, Ey = pyee.post.getField(data, sol, "Ey")
points0 = (Z[:, 0], R[0, :])
Ey = np.reshape(interpn(points0, Ey, points1, bounds_error=False, fill_value=0), np.shape(Z1))

# Compute power absorption
eps0 = 8.854187817e-12
omega = 2*np.pi*data['simulation']['freq']

Jp = np.zeros((Z1.shape[0], Z1.shape[1], 3))
Evect = np.zeros((Z1.shape[0], Z1.shape[1], 3))
Dvect = np.zeros((Z1.shape[0], Z1.shape[1], 3))

for iz in range(0, 2*data['geometry']['nz'] + 1):
    for ix in range(0, 2*data['geometry']['nx'] + 1):
        Evect[iz, ix, :] = np.array([Ez[iz, ix], Ex[iz, ix], Ey[iz, ix]], dtype=complex)
        Dvect[iz, ix, :] = np.inner(data['kappa'][iz, ix, :, :], Evect[iz, ix, :])
        Jp[iz, ix, :] = -1j*np.inner(data['kappa'][iz, ix, :, :] - np.eye(3), Evect[iz, ix, :])*omega*eps0

def divergence(field):
    "return the divergence of a n-D field"
    dz = np.max(Z1[:])/Z1.shape[0]
    dx = np.max(R1[:])/Z1.shape[1]

    dFZ = np.gradient(field[:,:,0], dz, dx)[0]
    dFX = np.gradient(field[:,:,1]*R1, dz, dx)[1]/R1
    dFY = 1j*data['simulation']['modes'][0]*field[:,:,2]/R1

    div = np.real(dFX + dFZ + dFY)
    return div

# pyee.post.plot(data, 'div(Jp) [A/m3]', Z = Z1, R = R1, scale = 'log', field=divergence(Jp))
# pyee.post.plot(data, 'div(Ja) [A/m3]', Z = Z1, R = R1, scale = 'log', field=divergence(data['current']))

Dref = np.mean(abs(Dvect))/np.sqrt(np.max(Z1[:])*np.max(R1[:]))
Eref = np.mean(abs(Evect))/np.sqrt(np.max(Z1[:])*np.max(R1[:]))
Dref = 1
Eref = 1

# pyee.post.plot(data, 'div(E)', Z = Z1, R = R1, field=abs(divergence(Evect))/Eref)
pyee.post.plot(data, 'div(D)', Z = Z1, R = R1, field=abs(divergence(Dvect)))
# pyee.post.plot(data, 'div(E)', Z = Z1, R = R1, scale = 'log', field=abs(divergence(Evect))/Eref)
# pyee.post.plot(data, 'div(D)', Z = Z1, R = R1, scale = 'log', field=abs(divergence(Dvect))/Dref)