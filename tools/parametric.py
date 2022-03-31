'''
Creation Date: Invalid date
Author: Pedro Jimenez, pejimene@ing.uc3m.es
##################################################
File: Untitled-1
Project: <<projectname>>
##################################################
Last Modified: Tuesday, April 27th 2021, 11:17:24 am
Modified By: Pedro Jimenez
##################################################
Code Description:
##################################################
Copyright (c) 2020 Plasma & Space Propulsion Team (EP2). Universidad Carlos III de Madrid.
'''
import run_pyee
import numpy as np

# # Helix
# def sim(h):

#     simulation = dict()
#     simulation['helix'] = float(h)
#     run_pyee.sim('HPT05M', simulation = simulation)

# # Mesh resolution
def nodes(case, nz, nx):

    geometry = dict()
    geometry['nz'] = nz
    geometry['nx'] = nx
    run_pyee.sim(case, geometry = geometry)



def cores(case, threads):
    general = dict()

    general['nthreads'] = threads
    run_pyee.sim(case, general = general)

