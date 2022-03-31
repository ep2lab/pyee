import scipy.sparse.linalg
import numpy as np
import pyee
import matplotlib.pyplot as plt
import h5py

file = 'helicon'

user_data = dict()

data = pyee.pre.prepareData('cases/' + file + '/', user_data)

freq = np.array([1, 5, 10, 50, 100, 238, 500, 1000, 5000, 10000])
cond = np.zeros(freq.shape)

for i in range(0, freq.size):

    data['simulation']['freq'] = freq[i]

    A, b = pyee.core.constructProblem(data, data['simulation']['modes'][0])

    norm_A = scipy.sparse.linalg.norm(A)
    norm_invA = scipy.sparse.linalg.norm(scipy.sparse.linalg.inv(A))
    cond[i] = norm_A*norm_invA

plt.plot(freq, np.log10(cond))

filename='condition.h5'
hf = h5py.File(filename, 'w')

hf.create_dataset('cond', data=cond)
hf.create_dataset('freq', data=freq)

hf.close()