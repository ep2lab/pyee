import sys
sys.path.append('../pyee/')
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as scp
import time
import pyee

t = time.time()

# Options
data = pyee.pre.prepareData()

A, b = pyee.core.constructProblem(data, data['simulation']['modes'][0])

solution = scp.spsolve(A, b)

elapsed = time.time() - t
print("Elapased time: %.4f seconds"% (elapsed))

# Analytical fields
freq = data['simulation']['freq']*1e6 # Recover dimensions
R, Z = np.meshgrid(np.linspace(0, data['geometry']['Lx']/(2*np.pi*freq/3e8), 2*data['geometry']['nx'] + 1),
                   np.linspace(0, data['geometry']['Lz']/(2*np.pi*freq/3e8), 2*data['geometry']['nz'] + 1))

p = plt.contourf(Z, R, abs(data['EA'][:, :, 0]), levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Ez')
plt.show()

p = plt.contourf(Z, R, abs(data['EA'][:, :, 1]), levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Er')
plt.show()

p = plt.contourf(Z, R, abs(data['EA'][:, :, 2]), levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Ephi')
plt.show()

# Current plots
p = plt.contourf(Z, R, data['current'][:, :, 0], levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Jz')
plt.show()

p = plt.contourf(Z, R, data['current'][:, :, 1], levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Jr')
plt.show()

p = plt.contourf(Z, R, abs(data['current'][:, :, 2]), levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Jphi')
plt.show()

# Error plots(Project to mesh)
Z, R, Ez = pyee.core.getField(data, solution, "Ez")
field = data['EA'][1:2*data['geometry']['nz'] + 1:2, 
                   0:2*data['geometry']['nz'] + 1:2, 0]
Ez = Ez[0:(np.shape(Ez)[0] - 1), :]
Z  = Z[0:(np.shape(Z)[0] - 1), :]
R  = R[0:(np.shape(R)[0] - 1), :]
p = plt.contourf(Z, R, abs(field - Ez)/np.max(abs(field))*100, levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Error Ez %')
plt.show()

Z, R, Ex = pyee.core.getField(data, solution, "Ex")
field = data['EA'][0:2*data['geometry']['nz'] + 1:2, 
                   1:2*data['geometry']['nz'] + 1:2, 1]
Ex = Ex[:, 0:(np.shape(Ex)[1] - 1)]
Z  = Z[:, 0:(np.shape(R)[1] - 1)]
R  = R[:, 0:(np.shape(R)[1] - 1)]
p = plt.contourf(Z, R, abs(field - Ex)/np.max(abs(field))*100, levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Error Er %')
plt.show()

Z, R, Ey = pyee.core.getField(data, solution, "Ey")
field = data['EA'][0:2*data['geometry']['nz'] + 1:2, 
                   0:2*data['geometry']['nz'] + 1:2, 2]
p = plt.contourf(Z, R, abs(field - Ey)/np.max(abs(field))*100, levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Error Ephi %')
plt.show()

# Plot numerical fields
Z, R, Ez = pyee.core.getField(data, solution, "Ez")
p = plt.contourf(Z, R, abs(Ez), levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Numerical Ez')
plt.show()

Z, R, Ex = pyee.core.getField(data, solution, "Ex")
p = plt.contourf(Z, R, abs(Ex), levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Numerical Er')
plt.show()

Z, R, Ey = pyee.core.getField(data, solution, "Ey")
p = plt.contourf(Z, R, abs(Ey), levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Numerical Ephi')
plt.show()