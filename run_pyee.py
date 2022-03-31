# """ Run PYEE Full Maxwell Finite Difference Solver"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import dtype
import pyee
import scipy.sparse.linalg as scp
import time
import pyee
import os
import sys
import pickle

# import condition_lib as cond

def sim(file='', mode = 'default', **user_data):

    print("Running Pyee")

    t = time.time()

    if not(file == ''):
        path = 'cases/' + file + '/'
        sys.path.append(path)
    else:
        path = './'

    # Options
    t1 = time.time()
    if mode == 'hyphen':
        data = pyee.hyphen.pre(path)
    else:
        data = pyee.pre.prepareData(path, mode, user_data)

    nz = data['geometry']['nz']
    nx = data['geometry']['nx']

    elapsed = time.time() - t1
    print("Preprocessor time: %.4f seconds"% (elapsed))
    data['general']['pre_time'] = elapsed

    if "fem" in data['general']:
        FEM = data['general']["fem"]
    else:
         FEM = 0

    # FEM = 0
    ################################################################################################################################
    """Solver"""

    if FEM:
        data['general']['assembly_time'] = 0.0
        t3 = time.time()

        import random
        nid = random.randrange(0,1e9)

        tempname = './temp' + str(nid) + '/'

        os.system('mkdir ' + tempname)
        os.system('cp ./' + data['general']['simdir'] + '/preproc.pkl ' + tempname + '/data.pkl')
        # os.system('rm FEM.log')

        nthreads = data['general']['nthreads']

        print('Running FEM solver with ' + str(nthreads) + ' threads')

        import subprocess   # Do not import fem.py in pyee __init__
        command = 'mpirun -np ' + str(nthreads) + ' python3 ./pyee/fem.py ' + tempname + ' ' + str(data['simulation']['modes'][0])
        with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True) as p, \
            open(data['general']['simdir'] + 'FEM.log', 'ab') as file:
            for line in p.stdout: # b'\n'-separated lines
                sys.stdout.buffer.write(line) # pass bytes as is
                file.write(line)

        elapsed = time.time() - t3
        print("Solve time: %.4f seconds"% (elapsed))
        data['general']['solve_time'] = elapsed

        # Reconstruct from parallel solution
        for core in range(0, nthreads):
            file     = open(tempname + '/sol'+str(core)+'.pkl',"rb")
            temp_sol = pickle.load(file)
            file.close()
            if core == 0:
                    fem_sol = temp_sol
            else:
                cond = np.logical_not(np.isnan(temp_sol['Ez']))
                for field_tag in ['Ez', 'Ex', 'Ey', 'Jz', 'Jx', 'Jy', 'Qa']:
                    fem_sol[field_tag][cond] = temp_sol[field_tag][cond]

        os.system('rm -rf ' + tempname)

        if data['geometry']['nz'] == data['mesh']['Z'].shape[0] - 1:  # New FEM meshes (not staggered):
            v = fem_sol
        else:                                                         # Quick fix for legacy FD staggered meshes
            v = np.zeros(6 * (nz + 1) * (nx + 1), dtype=complex)
            for iz in range(0,nz):
                for ix in range(0,nx):  
                    v[6 * iz + 6 * (nz + 1) * ix + 0] = fem_sol['Ez'][iz,ix]
                    v[6 * iz + 6 * (nz + 1) * ix + 1] = fem_sol['Ex'][iz,ix]
                    v[6 * iz + 6 * (nz + 1) * ix + 2] = fem_sol['Ey'][iz,ix]
    else:
        print('Running FD solver')
        t2 = time.time()
        A, b = pyee.assembly.maxwell.constructProblem(data, data['simulation']['modes'][0])
        elapsed = time.time() - t2
        print("Assembly time: %.4f seconds"% (elapsed))
        data['general']['assembly_time'] = elapsed

        t3 = time.time()
        omega = 2*np.pi*data['simulation']['freq']*1.0e6
        fact =  4*np.pi*1.0e-7*omega*(omega/3.0e8)**(-2)
        v = scp.spsolve(A, b*fact)
        elapsed = time.time() - t3
        print("Solve time: %.4f seconds"% (elapsed))
        data['general']['solve_time'] = elapsed

    # t4 = time.time()
    # # sigmaL = scp.svds(A, k=1, which = 'LM', return_singular_vectors=False)
    # # sigmaS = scp.svds(A, k=1, which = 'SM', return_singular_vectors=False)
    # kappa1 = cond.condition_linpack(A.shape[1], A)
    # elapsed = time.time() - t4
    # print("Condition time: %.4f seconds"% (elapsed))
    # print("{:e}".format(kappa1))
    elapsed = time.time() - t
    data['general']['time'] = elapsed
    ################################################################################################################################
    """ Postprocessing """
    t4 = time.time()

    if data['general']['export']:
        pyee.post.export(data, v)

    if mode == 'hyphen':
        pyee.hyphen.post(path, data, v)

    if os.path.isfile(path + 'post_script.py'):
        import post_script
        post_script.postprocess(data, v)
        print("Pyee Finished")
    
    elapsed = time.time() - t4
    print("Postprocessing Time: %.4f seconds"% (elapsed))

    elapsed = time.time() - t
    print("Total Computational Time: %.4f seconds"% (elapsed))
 
    return data, v

if __name__ == '__main__':
    # Map command line arguments to function arguments.

    sim(*sys.argv[1:])


# """Old check mode"""
#  if data['simulation']['check_mode'] == 'on':

#         Z, R, Ez = pyee.post.getField(data, v, "Ez")
#         analytic = data['EA'][1:2*data['geometry']['nz'] + 1:2, 
#                               0:2*data['geometry']['nz'] + 1:2, 0]*fact
#         pyee.post.plot(data, 'error Ez %', Z=Z, R=R, 
#                        field=(abs(analytic) - abs(Ez))/np.mean(abs(analytic))*100)

#         Z, R, Ex = pyee.post.getField(data, v, "Ex")
#         analytic = data['EA'][0:2*data['geometry']['nz'] + 1:2, 
#                               1:2*data['geometry']['nz'] + 1:2, 1]*fact
#         pyee.post.plot(data, 'error Ex %', Z=Z, R=R, 
#                        field=(abs(analytic) - abs(Ex))/np.mean(abs(analytic))*100)

#         Z, R, Ey = pyee.post.getField(data, v, "Ey")
#         analytic = data['EA'][0:2*data['geometry']['nz'] + 1:2, 
#                               0:2*data['geometry']['nz'] + 1:2, 2]*fact
#         pyee.post.plot(data, 'error Ey %', Z=Z, R=R, 
#                        field=(abs(analytic) - abs(Ey))/np.mean(abs(analytic))*100)

#         plt.show()