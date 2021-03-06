# PYEE

Pyee solves the **frequency domain Maxwell equations** in a planar or axisymmetric domain using either a mixed element FEM formulation or a FD Yee discretization scheme. The new FEM solver is implemented through the parallel FEniCS/dolfinx libary and modern direct linear PETSc solvers. Both the assembly and solver steps are parallelized

The medium is characterized by a general complex dielectric tensor. The current version supports rectangular and uniform meshes with arbitrary integer **azimuthal mode expansions** in the axissymetric mode and **metallic walls (PEC)** boundary conditions. The code was developed for simulations of the **electromagnetic wave propagation through magnetized cold plasmas** and thus includes additional functions for this analysis

## Installation

***NOTE: Remember to export your $PYEE enviroment variable*** For example if the folder is located at the user folder and bash is used as terminal:
    
    echo 'export PYEE=$HOME/pyee' >> $HOME/.bashrc

The only current installation option is using containers until DOLFINx hits a stable version
  
###  Using docker containers
The FEM module contains a considerable number of dependencies and precompiled libraries. Containers are a great option for an easy installation of the code and reproductability across different machines and architectures. You will need to install the Docker engine https://docs.docker.com/engine/install/ (needs sudo permissions) and make sure to follow the post-installation steps to 'Manage Docker as a non-root user' and test docker

    sudo groupadd docker
    sudo usermod -aG docker $USER
    newgrp docker
    docker run hello-world
    
Future versions will include the experimental Docker rootless installation https://docs.docker.com/engine/security/rootless/
Download the public ``ep2lab\pyee`` image from Docker Hub and run the installation script (this will only match your user and group IDs inside the image so you can modify outputs files with user permissions later:
  
    docker pull ep2lab/pyee:latest
    cd docker
    ./install.sh

Now you have 2 options. You can forget about the container and run the specialized script that will mount the filesystem and manage the container under the hood (see next section). Or you can work inside a container:

      docker run -ti ep2lab/pyee:latest
      
This is often useful for development and debugging or necessary if you are running certain OS (for example windows). You can get familiar with many Docker commands taking a look at the cheatsheet included in the ``docs/`` folder. Find more information about building images at https://hub.docker.com/r/ep2lab/pyee 

***NOTE: Using the script ``pyee_docker.sh`` will mount the source code from ``$PYEE`` folder. Working inside the container will execute the code delivered at the time the image was pushed to Docker Hub (unless you mount a file system manually)***

## Simulation Setup, Pre and Postprocessing Options

All setup parameters and postprocessing calls are provided in the case folder ``cases/<<case_name>>``

The simulation configuration is given in the ``simrc.yml`` file. This also provides solver (FEM or FD) and controls the main simulation parameters

The current density and dielectric tensor can be given by an optional python script (named ``config_sim.py``) inside each case folder or as simple functions in ``simrc.yml``. The preprocessor can compute the dielectric tensor from plasma density and collisionality profiles + the background magnetic field topology. You can pass additional user defined arguments in ``simrc.yml`` that can be access in ``config_sim.py`` through the internal ``data`` dictionary! This is very useful for quick configuration changes and parametric analysis

Several postprocessing functions, such as plotting and diagnostic utilities, are provided in the ``post.py`` module. These functions can be called in  ``post_script.py`` included in each case folder ``cases/<<case_name>>``. The ``post_script.py`` is optional and automatically detected by the code

There is also an option in ``simrc.yml`` to export the electric and magnetic fields and other relevant data such as the power profile and the computation time and outputs a file in the .h5 format (default ``result.h5``) without using any ``post_script.py``

## Run simulation case

Examples are provided in the ``cases/`` folder. To run via terminal do:

    ./pyee.sh <<case_name>> <<mode>>
  
To run the code via a preinstalled container simply use:
  
    ./pyee_docker.sh <<case_name>> <<mode>>
  
The available modes keywords are:

    default    -> runs preprocessing and solver 
    
    solve_pre  -> loads 'preproc.pkl' and solves. This overwrittes all the setups in 'simrc.yml' except the 'general' section
    
    hyphen     -> runs HYPHEN wrapper pre and post functions

Inside Python or IPython:

    data, solution = run_pyee.sim('<<case_name>>', mode = <<mode>>, **user_data) 
  
This also outputs the preprocessing and solution data structures. ``user_data`` is an optional dictionary that overwrittes ``simrc.yml`` fields

## Additional Tools
  
Inside the ``tools/`` folder, there is a whole ``matlabPost`` MATLAB set of functions meant for postprocessing and plotting but also capable of calling the PYEE module for execution **(experimental, not documented)**

The ``antenna.py`` function computes the current density for typical antenna setups. Currently only the helix family of antennas is available. The antenna configuration is given in the dictionary ``antenna_data`` and the azimuthal mode for the simulation must be provided:

    import tools.antenna as antenna
    antenna.helix(data, antenna_data, mode)
 
The code can also run parametric studies modifying some configuration parameters. Look at ``run_parametric.py`` and ``tools/parametric.py`` to learn how
