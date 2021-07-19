Getting Started
===============

########################
Introduction 
########################

Gaussian accelerated molecular dynamics (GaMD) is a well-established enhanced sampling method for molecular dynamics (MD) simulations that effectively samples the potential energy landscape of the system by adding a boost potential, which smoothens the surface and lowers energy barriers between states. However, GaMD is unable to give time-dependent properties like kinetics directly. On the other hand, the weighted ensemble (WE) method can efficiently sample transitions between states with its many weighted trajectories, which directly give us rates and pathways. However, the WE method's performance (i.e., convergence and efficiency) depends heavily on its initial conditions or initial sampling of the potential energy landscape. Hence, we have developed a hybrid method that combines the two methods. GaMD is first used to sample the potential energy landscape of the system. Then the WE method is used to sample further the potential energy landscape and kinetic properties of interest. We show that the hybrid method can sample both thermodynamic and kinetic properties more accurately and quickly compared to using one method by itself. gamd_we aims to create starting structures for WE simulations.

########################
Software Requirements
########################

* Amber 18
* Cuda

########################
Installation and Setup Instructions
########################

* Make sure `anaconda3 <https://www.anaconda.com/>`_ is installed on the local machine. 
* Go to the `download <https://www.anaconda.com/products/individual>`_  page of anaconda3 and install the latest version of anaconda3. 
* Create a new conda environment with python = 3.8 and install the package with the following commands in the terminal: 

.. code-block:: python

    conda create -n gamdwe python=3.6 # Create a new conda environment

.. code-block:: python

    conda activate gamdwe # Activate the conda environment

.. code-block:: python

    conda install openforcefield # Install openforcefield

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

    conda install -c conda-forge jupyterlab # Install jupyterlab

.. code-block:: python

    conda install git # Install git

* Clone the *gamd_we* repository :

.. code-block:: python

    git clone https://github.com/anandojha/gamd_we.git

########################
Gaussian Accelerated Molecular Dynamics
########################


