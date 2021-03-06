Getting Started
===============

################################################################################################
Introduction 
################################################################################################

Gaussian accelerated molecular dynamics (GaMD) is a well-established enhanced sampling method for molecular dynamics (MD) simulations that effectively samples the potential energy landscape of the system by adding a boost potential, which smoothens the surface and lowers energy barriers between states. However, GaMD is unable to give time-dependent properties like kinetics directly. On the other hand, the weighted ensemble (WE) method can efficiently sample transitions between states with its many weighted trajectories, which directly give us rates and pathways. However, the WE method's performance (i.e., convergence and efficiency) depends heavily on its initial conditions or initial sampling of the potential energy landscape. Hence, we have developed a hybrid method that combines the two methods. GaMD is first used to sample the potential energy landscape of the system. Then the WE method is used to sample further the potential energy landscape and kinetic properties of interest. We show that the hybrid method can sample both thermodynamic and kinetic properties more accurately and quickly compared to using one method by itself. 

gamd_we aims to create starting structures for WE simulations. A series of functions in the module downloads the required PDB (alanine dipeptide and chignolin systems), prepare the system in OpenMM (equilibration followed by solvation) and run GaMD simulations in Amber followed by reweighing, analysis and generating starting structures for Weighted Ensemble Molecular Dynamics. We need CUDA enabled versions of OpenMM and Amber to run the scripts (or necessary changes can be made in these scripts, e.g., for the platform keyword in the OpenMM codes and pmemd keywords in Amber) to run on non-CUDA machines.

################################################################################################
Software Requirements
################################################################################################

* Amber 20
* Cuda

################################################################################################
Installation and Setup Instructions
################################################################################################

* Make sure `anaconda3 <https://www.anaconda.com/>`_ is installed on the local machine. 
* Go to the `download <https://www.anaconda.com/products/individual>`_  page of anaconda3 and install the latest version of anaconda3. 
* Create a new conda environment with python = 3.6 and install the package with the following commands in the terminal: 

.. code-block:: python

    conda create -n gamdwe python=3.6 # Create a new conda environment

.. code-block:: python

    conda activate gamdwe # Activate the conda environment

.. code-block:: python

    conda install -c conda-forge curl # Install curl

.. code-block:: python

    conda install -c conda-forge matplotlib # Install matplotlib

.. code-block:: python

    conda install -c conda-forge openmm  # Install openmm

.. code-block:: python

    conda install -c conda-forge seaborn # Install seaborn

.. code-block:: python

    conda install -c conda-forge pandas # Install pandas

.. code-block:: python

    conda install -c conda-forge mdtraj # Install mdtraj

.. code-block:: python

    conda install -c conda-forge openmm # Install openmm

.. code-block:: python

    conda install git # Install git

* Clone the *gamd_we* repository :

.. code-block:: python

    git clone https://github.com/anandojha/gamd_we.git

################################################################################################
Gaussian Accelerated Molecular Dynamics and Amber
################################################################################################

Gaussian Accelerated Molecular Dynamics (GaMD) has been implemented in pmemd, both the serial and parallel versions on CPU (pmemd and pmemd.MPI) and  GPU (pmemd.cuda and pmemd.cuda.MPI). Similar to AMD, GaMD provides options to add only the total potential boost (igamd=1), only dihedral potential boost (igamd=2), or the dual potential boost (igamd=3). The dual-boost simulation generally provides higher acceleration than the other two types of simulations for enhanced sampling. Simulation parameters comprise of the threshold energy values and the effective harmonic force constants,  k\ :sub:`0P`\   and   k\ :sub:`0D`\  for the total and dihedral potential boost, respectively. All information generated by GaMD necessary for reweighing is stored at each step into a vector which is flushed to a log file (gamd.log by default) every time the coordinates are written to disk, i.e., every two steps. The following variables specify additional parameters.

1. imin  -  Flag to run minimization 

.. code-block:: none

    = 0 (default) run molecular dynamics without any minimization

.. code-block:: none

    = 1 perform an energy minimization

.. code-block:: none

    = 5 read in a trajectory for analysis


2. irest  - Flag to restart simulation 

.. code-block:: none 

    = 0 (default) do not restart the simulation and run as a new simulation

.. code-block:: none 

    = 1 restart the simulation reading coordinates and velocities from a previously saved restart file

3. ntx  -  Option to read initial coordinates, velocities and box size from the inpcrd file 

.. code-block:: none 

    = 1 (default) coordinates but no velocities will be read from either a NetCDF or a formatted (ASCII) coordinate

.. code-block:: none 

    = 5 coordinates and velocities will be read from either a NetCDF or a formatted (ASCII) coordinate. Box information will be read if ntb > 0 and velocity information will only be used if irest = 1

4. nstlim  - Number of MD steps to be performed

5. dt  - Timesteps (in picoseconds) 

6. ntc  -  Flag for SHAKE to perform bond length constraints 

.. code-block:: none 

    = 1 SHAKE is not performed (default) 

.. code-block:: none 

    = 2 bonds involving hydrogen are constrained 

.. code-block:: none 

    = 3 all bonds are constrained 


7. ntf - To employ TIP3P, set ntf = NTC = 2

8. tol - Tolerance for convergence

9. wrap   -  If iwrap = 1, coordinates written to the restart and trajectory files will be wrapped into a primary box. This means that for each molecule, its periodic image is closest to the middle of primary box (with x coordinates between 0 and a, y coordinates between 0 and b, and z coordinates between 0 and c) will be the one written to the output file. This often makes resulting structures look better visually, but has no effect on the energy or forces. Performing such wrapping, however, can mess up diffusion and other calculations. If iwrap = 0, no wrapping will be performed, in which case it is typical to use cpptraj as a post-processing program to translate molecules back to the primary box. For very long runs, setting iwrap = 1 may be required to keep the coordinate output from overflowing the trajectory and restart file formats, especially if trajectories are written in ASCII format instead of NetCDF.


10. ntb - This variable controls whether or not periodic boundaries are imposed on the system during calculation of non-bonded interactions. Bonds spanning periodic boundaries are not yet supported. There is no longer any need to set this variable, since it can be determined from igb and ntp parameters. The ???proper??? default for ntb is chosen (ntb = 0 when igb > 0, ntb = 2 when ntp > 0, and ntb = 1 otherwise). This behavior can be overridden by supplying an explicit value, although this is discouraged to prevent errors. 

.. code-block:: none 

    = 0 no periodicity is applied and PME is off (default when igb > 0)

.. code-block:: none 

    = 1 constant volume (default when igb and NTP are both 0, which are their defaults)

.. code-block:: none 

    = 2 constant pressure (default when ntp > 0) 


11. cut - Non bonded cut-off distance

12. temp0 - Reference temperature at which the system is to be kept, if not > 0. For temperatures above 300 K, the step size should be reduced since increased distance traveled between evaluations can lead to SHAKE and other problems. (Default 300 K)

13. gamma_ln - Collision frequency for Langevin dynamics, in 1/ps. Values in the range 2/ps - 5/ps often give acceptable temperature control, while allowing transitions to take place. Values near 50/ps correspond to the collision frequency for liquid water, and maybe useful if rough physical time scales for motion are desired.

14. ntpr - Every ntpr steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo". "mdinfo" is closed and reopened each time, so it always contains the most recent energy and temperature. (Default 50)

15. ntwx  - Every ntwx steps, the coordinates will be written to the mdcrd file. If ntwx = 0, no coordinate trajectory file will be written. (Default = 0)

16. ntwr  - Every two steps during dynamics, the ???restart??? file will be written, ensuring that recovery from a crash will not be so painful. No matter what the value of ntwr, a restart file will be written at the end of the run, i.e., after nstlim steps (for dynamics) or maxcyc steps (for minimization). If ntwr < 0, a unique copy of the file, restrt_nstep, is written every abs(ntwr) steps. This option is useful if, for example, one wants to run free energy perturbations from multiple starting points or save a series of restart files for minimization. Default = nstlim.

17. ntxo - Format of the final coordinates, velocities, and box size (if constant volume or pressure run) written to file "restart" 

.. code-block:: none 

    = 1 Formatted (ASCII)

.. code-block:: none 
  
    = 2 (default) NetCDF file


18. ioutfm - Format of coordinate and velocity trajectory files (mdcrd, mdvel and inptraj). Binary output is in NetCDF trajectory format 

.. code-block:: none 

    = 0 Formatted ASCII trajectory

.. code-block:: none 

    = 1 (default) Binary NetCDF trajectory 


19. in - Random seed


20. ntwprt  - Number of atoms to include in trajectory files (mdcrd and mdvel). This flag can be used to decrease the size of the these files, by including only the first part of the system, which is usually of greater interest (for instance, one might include only the solute and not the solvent) 

.. code-block:: none 

    = 0 (default) Include all atoms of the system when writing trajectories.


21. igamd  - Flag to apply boost potential 

.. code-block:: none  

    = 0 (default) no boost is applied 

.. code-block:: none  

    = 1 boost on the total potential energy only

.. code-block:: none 

    = 2 boost on the dihedral energy only

.. code-block:: none 

    = 3 dual boost on the both, dihedral and total potential energy 


22. iE - Flag to set the threshold energy E 

.. code-block:: none 

    = 1 (default) set the threshold energy to the lower bound E = Vmax

.. code-block:: none 

    = 2 set the threshold energy to the upper bound E = Vmax + 1/k 

23. irest_gamd - Flag to restart GaMD simulation 

.. code-block:: none 

    = 0 (default) new simulation. A file "gamd-restart.dat" that stores the maximum, minimum, average and standard deviation of the dihedral and/or total potential energies (depending on the igamd flag) will be saved automatically after GaMD equilibration stage

.. code-block:: none

    = 1 restart simulation (ntcmd is set to 0 in this case). The "gamd-restart.dat" file will be read for restart 


24. ntcmd - Number of initial conventional molecular dynamics simulation steps. Potential energies are collected between ntcmdprep and ntcmd to calculate their maximum, minimum, average and standard deviation (V\ :sub:`max`\, V\ :sub:`min`\, V\ :sub:`avg`\, ??\ :sub:`V`\). The default is 1,000,000 for a simulation with 2 fs timestep.


25. nteb - Number of biasing molecular dynamics simulation steps. Potential statistics (V\ :sub:`max`\, V\ :sub:`min`\, V\ :sub:`avg`\,  ??\ :sub:`V`\) are updated between the ntebprep and nteb steps and used to calculate GaMD acceleration parameters, particularly E and  k\ :sub:`0`\. The default is 1,000,000 for a simulation with 2 fs timestep. A greater value may be needed to ensure that the potential statistics and GaMD acceleration parameters level off before running production simulation between the nteb and nstlim (total simulation length) steps. Moreover, nteb can be set to nstlim, by which the potential statistics and GaMD acceleration parameters are updated adaptively throughout the simulation. This is some cases provides more appropriate acceleration.

26. ntave - Number of simulation steps used to calculate the average and standard deviation of potential energies. The default is set to 50,000 for GaMD simulations. It is recommended to be updated as about 4 times of the total number of atoms in the system. Note that ntcmd and nteb need to be multiples of ntave.

27. ntcmdprep -  Number of preparation conventional molecular dynamics steps. This is used for system equilibration and the potential energies are not collected for calculating their statistics. The default is 200,000 for a simulation with 2 fs timestep.

28. ntebprep - Number of preparation biasing molecular dynamics simulation steps. This is used for system equilibration after adding the boost potential and the potential statistics  (V\ :sub:`max`\, V\ :sub:`min`\, V\ :sub:`avg`\,  ??\ :sub:`V`\) are not updated during these steps. The default is 200,000 for a simulation with 2 fs timestep.

29. sigm0D - Upper limit of the standard deviation of the dihedral potential boost that allows for accurate reweighting if igamd is set to 2 or 3. The default is 6.0 (unit: kcal/mol).
30. sigma0P - Upper limit of the standard deviation of the total potential boost that allows for accurate reweighting if igamd is set to 1 or 3. The default is 6.0 (unit: kcal/mol).
