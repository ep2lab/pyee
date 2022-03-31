import sys
sys.path.append('')
import pyee
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as scp
import time

t = time.time()

data = pyee.pre.prepareData()
A, b = pyee.core.constructProblem(data, 1)

solution = scp.spsolve(A, b)

elapsed = time.time() - t
print("Elapased time: %.4f seconds"% (elapsed))

freq = data['simulation']['freq']*1e6 # Recover dimensions
R, Z = np.meshgrid(np.linspace(0, data['geometry']['Lx']/(2*np.pi*freq/3e8), 2*data['geometry']['nx'] + 1),
                   np.linspace(0, data['geometry']['Lz']/(2*np.pi*freq/3e8), 2*data['geometry']['nz'] + 1))
# Current plots
p = plt.contourf(Z, R, abs(data['current'][:, :, 0]), levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Jz')
plt.show()

p = plt.contourf(Z, R, abs(data['current'][:, :, 2]), levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title('Jz')
plt.show()

# Plasma plots
p = plt.contourf(Z, R, abs(data['plasma'][:, :, 0]), levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title(r'$\omega_{ce}$')
plt.show()

p = plt.contourf(Z, R, abs(data['plasma'][:, :, 2]), levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title(r'$\omega_{pe}$')
plt.show()

p = plt.contourf(Z, R, abs(data['plasma'][:, :, 6]), levels = 100)
plt.colorbar(p)
for c in p.collections:
    c.set_edgecolor("face")
plt.title(r'$\theta$')
plt.show()

# p = plt.streamplot(Z[:, 0], R[0, :], abs(np.cos(data['plasma'][:, :, 6])*data['plasma'][:, :, 2]**2), 
#                     abs(np.sin(data['plasma'][:, :, 6])*data['plasma'][:, :, 2]**2), density = 1)
# plt.title('Magenetic field lines')

# Field plots
Z, R, Fz = pyee.core.getField(data, solution, "Ez")
fig, ax = plt.subplots()
cs   = ax.contourf(Z, R, abs(Fz), levels = 100)
cbar = fig.colorbar(cs)
for c in cs.collections:
    c.set_edgecolor("face")
plt.title('Ez')
plt.show()

Z, R, Fx = pyee.core.getField(data, solution, "Ex")
fig, ax = plt.subplots()
cs   = ax.contourf(Z, R, abs(Fx), levels = 100)
cbar = fig.colorbar(cs)
for c in cs.collections:
    c.set_edgecolor("face")
plt.title('Ex')
plt.show()

Z, R, Fy = pyee.core.getField(data, solution, "Ey")
fig, ax = plt.subplots()
cs   = ax.contourf(Z, R, abs(Fy), levels = 100)
cbar = fig.colorbar(cs)
for c in cs.collections:
    c.set_edgecolor("face")
plt.title('Ey')
plt.show()

# LOG PLOTS

Z, R, Fz = pyee.core.getField(data, solution, "Ez")
fig, ax = plt.subplots()
cs   = ax.contourf(Z, R, np.log(abs(Fz)), levels = 100)
cbar = fig.colorbar(cs)
for c in cs.collections:
    c.set_edgecolor("face")
plt.title('Log Ez')
plt.show()

Z, R, Fx = pyee.core.getField(data, solution, "Ex")
fig, ax = plt.subplots()
cs   = ax.contourf(Z, R, np.log(abs(Fx)), levels = 100)
cbar = fig.colorbar(cs)
for c in cs.collections:
    c.set_edgecolor("face")
plt.title('Log Ex')
plt.show()

Z, R, Fy = pyee.core.getField(data, solution, "Ey")
fig, ax = plt.subplots()
cs   = ax.contourf(Z, R, np.log(abs(Fy)), levels = 100)
cbar = fig.colorbar(cs)
for c in cs.collections:
    c.set_edgecolor("face")
plt.title('Log Ey')
plt.show()

# TODO: Properly project to a common grid
# fig, ax = plt.subplots()
# cs   = ax.contourf(Z, R, (abs(Fx)**2 + abs(Fy)**2 + abs(Fz)**2)**0.5, levels = 30)
# cbar = fig.colorbar(cs)
# for c in cs.collections:
#     c.set_edgecolor("face")
# plt.title('Enorm')
# plt.show()