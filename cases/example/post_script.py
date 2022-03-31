import numpy as np
import pyee
import matplotlib.pyplot as plt
import h5py

def postprocess(data, solution):

    R1 = data['mesh']['R']
    Z1 = data['mesh']['Z']

    # pyee.post.plot(data, 'Jz', Z=Z1, R=R1, field=abs(data['current'][:, :, 0]))
    # pyee.post.plot(data, 'Jr', Z=Z1, R=R1, field=abs(data['current'][:, :, 1]))
    # pyee.post.plot(data, 'Jphi', Z=Z1, R=R1, field=abs(data['current'][:, :, 2]))

    # pyee.post.plot(data, 'omega_ce', Z=Z1, R=R1, field=data['plasma'][:, :, 0])
    # pyee.post.plot(data, 'omega_pe', Z=Z1, R=R1, field=data['plasma'][:, :, 2])
    # pyee.post.plot(data, 'theta', Z=Z1, R=R1, field=data['plasma'][:, :, 6])

    # p = plt.streamplot(Z[:, 0], R[0, :], abs(np.cos(data['plasma'][:, :, 6])*data['plasma'][:, :, 2]**2), 
    #                     abs(np.sin(data['plasma'][:, :, 6])*data['plasma'][:, :, 2]**2), density = 1)
    # plt.title('Magenetic field lines')

    # # # Plot fields
    plt.figure(figsize=(8, 6))
    pyee.post.plot(data, 'Ez', field_tag='Ez', solution=solution, plot_type='magnitude', saveFlag='on', levels=500)
    plt.figure(figsize=(8, 6))
    pyee.post.plot(data, 'Ex', field_tag='Ex', solution=solution, plot_type='magnitude', saveFlag='on', levels=500)
    plt.figure(figsize=(8, 6))
    pyee.post.plot(data, 'Ey', field_tag='Ey', solution=solution, plot_type='magnitude', saveFlag='on', levels=500)
    plt.figure(figsize=(8, 6))
    pyee.post.plot(data, 'Ez', field_tag='Ez', solution=solution, plot_type='phase', saveFlag='on', levels=500)
    plt.figure(figsize=(8, 6))
    pyee.post.plot(data, 'Ex', field_tag='Ex', solution=solution, plot_type='phase', saveFlag='on', levels=500)
    plt.figure(figsize=(8, 6))
    pyee.post.plot(data, 'Ey', field_tag='Ey', solution=solution, plot_type='phase', saveFlag='on', levels=500)
    # pyee.post.plot(data, 'Ez', field_tag='Ez', solution=solution, plot_type='phase', saveFlag='off')
    # pyee.post.plot(data, 'Ex', field_tag='Ex', solution=solution, plot_type='phase', saveFlag='off')
    # pyee.post.plot(data, 'Ey', field_tag='Ey', solution=solution, plot_type='phase', saveFlag='off')

    # pyee.post.plot(data, 'Ez', field_tag='Ez', solution=solution, scale='log', saveFlag='off')
    # pyee.post.plot(data, 'Ex', field_tag='Ex', solution=solution, scale='log', saveFlag='off')
    # pyee.post.plot(data, 'Ey', field_tag='Ey', solution=solution, scale='log', saveFlag='off')

    # # Compute Power
    # P = pyee.post.power_abs(data, solution)
    # P[P<1] = np.nan
    # # pyee.post.plot(data, 'Power', Z=Z1, R=R1, field=P)
    # pyee.post.plot(data, 'Power', Z=Z1, R=R1, field=P ,saveFlag='off')
    # pyee.post.plot(data, 'Power', Z=Z1, R=R1, field=P, scale='log', saveFlag='off')

    # Save data
    # pyee.post.export(data,solution, filename = 'results.h5')

    # Plot propagation regions
    # F = pyee.post.cma(Z1, R1, data['plasma'][:, :, 0], data['plasma'][:, :, 2], 0) #1/(1836*40)
    
    plt.show()
    return