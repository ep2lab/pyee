=======
PYEE
=======

Pyee solves the **frequency domain Maxwell equations** in a planar or axisymmetric domain using the Yee discretization scheme. 
The medium is characterized by a general complex dielectric tensor. The current version supports rectangular and uniform meshes with arbitrary integer **azimuthal mode expansions** in the axissymetric mode and **metallic walls (PEC)** boundary conditions. The code was developed for simulations of the **electromagnetic wave propagation through magnetized cold plasmas** and thus includes additional functions for this analysis. 

Simulation Setup
=======
All setup parameters and postprocessing calls are provided in the case folder ``cases/<<case_name>>``

The main simulation configuration is given in the ``simrc.yml`` file. The current density and dielectric tensor can be given by a python script (default name ``config_sim.py``) inside each case folder or as simple functions in ``simrc.yml``. The preprocessor can compute the dielectric tensor from plasma density and collisionality profiles + the background magnetic field topology (see Helicon case) through the ``pre.py kappa`` function.


Run simulation case
=======
Examples are provided in the ``cases/`` folder. To run via terminal do ::

  python run_pyee <<case_name>> 

Or inside Python or IPython::

  run_pyee.sim('<<case_name>>') 
  

Run simulation core
=======
For integration with external tools and other application a minimal wrapper is provided in the ``run_core.py`` file. This is the basic assembly and solver functionality with minimal preprocessing and postprocessing. The output are the electric fields and power in the PYEE mesh. The ``data`` configuration dictionary is provided as input indicating the azimuthal mode for axissymetric simulations::

  import run_core
  Z, R, Ez, Ex, Ey, P = run_core.sim(data, mode)
  
Postprocessing
=======
Several postprocessing functions, such as plotting and diagnostic utilities, are provided in the ``post.py`` module. These functions will be called in  ``post_script.py`` included in each case folder ``cases/<<case_name>>``.

The function ``save_hdf5`` gets the electric and magnetic fields and other relevant data such as the power profile and the computation time and outputs a file in the .h5 format (default ``result.h5``)`.

Additional Tools
=======
Inside the ``tools/`` folder, there is a whole ``matlabPost`` MATLAB set of functions meant for postprocessing and plotting but also capable of calling the PYEE module for execution **(experimental, not documented)**. 

The ``antenna.py`` function computes the current density for typical antenna setups. Currently only the helix family of antennas is available. The antenna configuration is given in the dictionary ``antenna_data`` and the azimuthal mode for the simulation must be provided::

  import antenna
  antenna.helix(data, antenna_data, mode)
  
**Notice that you might have to add** ``tools`` **to the path or copy these functions to the** ``cases/<<case_name>>`` **folder.**
