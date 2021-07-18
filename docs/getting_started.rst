Getting Started
===============

This page details how to get started with QMMMReBind. 

########################
Introduction 
########################


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

    conda install -c ambermd pytraj # Install pytraj

.. code-block:: python

    conda install -c omnia ambertools # Install ambertools

.. code-block:: python

    conda install git # Install git

* Clone the *gamd_we* repository :

.. code-block:: python

    git clone https://github.com/anandojha/gamd_we.git

########################
Notes on Gaussian Accelerated Molecular Dynamics
########################


