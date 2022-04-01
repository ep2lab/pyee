#!/bin/bash
# Inputs:    1) Simulation name (cases folder)
#            2) Simulation mode. 'default'    runs preprocessing and solver, 
#                                'solve_pre'  (loads preproc.pkl and solves), 
#                                'hyphen'     runs HYPHEN wrapper pre and post functions
# A valid image of FEniCS/dolfinx shoud be installed. Visit https://github.com/FEniCS/dolfinx

# REMEMBER TO EXPORT THE PYEE ENVIROMENT VARIABLE!!!! (example where pyee is in user folder):
# echo 'export PYEE=$HOME/pyee' >> $HOME/.bashrc

################################## INSTRUCTIONS FOR CONTAINERS #######################################
# 0) Install docker (sudo) https://docs.docker.com/engine/install/, "DON'T forget Manage Docker as a non-root user"
# 1) We recommend to install the standard dolfinx image (follow the instructions at github)
# 2) Then start interactive container 'docker run -ti --rm --name pyee_cont ep2lab/pyee:latest' (or the name of the image) 
# 3) Inside the container use pip3 to install additional dependencies (scipy, yaml, h5py, matplotlib)
# 4) In another terminal OUTSIDE container get your user and group ids, 'id -u' and 'id -g', write down this numbers
# 5) INSIDE image run this (Recommended so created files have user permissions):
        # addgroup --gid $GROUP_ID user
        # adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user
# 6) Now save the running container as a new image for example ep2lab/pyee:local 
        # docker commit pyee_cont ep2lab/pyee:local

# Now to RUN pyee:
# default mode:
# docker run --rm -v $PYEE:/root -u $(id -u):$(id -g) ep2lab/pyee:local /bin/bash pyee.sh example
# Or in hyphen mode:
# docker run --rm -v $PYEE:/root -u $(id -u):$(id -g) ep2lab/pyee:local /bin/bash pyee.sh HIP01_w hyphen
# Or solve  pre-processed simulation:
# docker run --rm -v $PYEE:/root -u $(id -u):$(id -g) ep2lab/pyee:local /bin/bash pyee.sh HIP01_w solve_pre

########################### Additional TIPS for containers (cheatshet at docs/) ####################:
# List available images                           -> docker images ls
# List running/all containers                     -> docker ps' or 'docker ps- a
# Eliminate stopped containers                    -> docker container prune
# Run an interactive container and delete at exit -> docker run -ti --rm
# Mount file systems                              -> docker run -v origin:dest
# If there is some kind of error inside the wave code you might need to prune containers !!!!!!!!!!!!!!!
#####################################################################################################

docker run --rm -v $PYEE:/root -u $(id -u):$(id -g) ep2lab/pyee:local /bin/bash pyee.sh $1 $2