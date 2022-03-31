#!/bin/bash
# Inputs:    1) Simulation name (cases folder)
#            2) Simulation mode. 'default'    runs preprocessing and solver, 
#                                'solve_pre'  (loads preproc.pkl and solves), 
#                                'hyphen'     runs HYPHEN wrapper pre and post functions
# A valid image of FEniCS/dolfinx shoud be installed. Visit https://github.com/FEniCS/dolfinx

# REMEMBER TO EXPORT THE PYEE ENVIROMENT VARIABLE!!!! (example where pyee is in user folder):
# echo 'export PYEE=$HOME/pyee' >> $HOME/.bashrc

source /usr/local/bin/dolfinx-complex-mode   # Updates enviroment variables
python3 -u run_pyee.py $1 $2                 # -u to prevent python buffering the stdout