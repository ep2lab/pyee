import numpy as np
import pyee
import run_pyee
import os
import shutil
import matplotlib.pyplot as plt

 
def copyDirectory(src, dest):
    try:
        shutil.copytree(src, dest)
    # Directories are the same
    except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        print('Directory not copied. Error: %s' % e)

points = 10

nz_ini = 107
nr_ini = 57

nz_end = 807
nr_end = 407

fz = (nz_end/nz_ini)**(1/points)
fr = (nr_end/nr_ini)**(1/points)

plt.figure()

user_data = dict()
user_data['geometry'] = dict()

# case = 'HPT03_F0'
for case in ('HPT03_F0','HPT03_F1','HPT03_F10','HPT03_F100','HPT03_F1000'):
    print(case)

    nnodes  = np.zeros(points + 1)
    Check  = np.zeros(points + 1)

    for i in range(0, points + 1):
        print(i)
        nz = round(nz_ini*fz**i)
        nr = round(nr_ini*fr**i)

        user_data['geometry']['nz'] = nz
        user_data['geometry']['nx'] = nr

        data, solution = run_pyee.sim(case, **user_data)

        if not(os.path.isdir(data['general']['simdir'] + '/convg_' + case)):
            os.mkdir(data['general']['simdir'] + '/convg_' + case)

        Z1, R1, Ez, Ex, Ey, P, time = pyee.post.save_hdf5(data,solution, filename = 'convg_' + case + '/results' + str(i) +'.h5')
        
        nnodes[i] = P.size
        Enorm = np.abs(Ex**2 + Ey**2 + Ez**2)
        Enorm[np.isnan(Enorm)] = 0

        Enorm[np.isnan(Enorm)] = 0

        Check[i] = np.max(Enorm)

    # plt.semilogy(nnodes, np.abs(Check[points] - Check)/Check[points], label = case)
    # copyDirectory('convg_' + case,'../temp3/'+'convg_' + case)

# plt.legend()
# plt.show()