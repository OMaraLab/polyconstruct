Getting Started
===============

This page details how to get started with PolyConstruct.


Setup for PolyConstruct:

From your home directory, install PolyConstruct from Git:

.. code-block:: python

    git clone https://github.com/OMaraLab/polyconstruct.git

Then navigate to polyconstruct:

.. code-block:: python

    cd polyconstruct

To setup polyconstruct, run: 

.. code-block:: python

    conda create --name polyconstruct-env
    conda activate polyconstruct-env

    pip install -r requirements.txt

    # polytop requires python 3.10
    # you may need to install python=3.10 before you are able to install requirements (as per above), depending on your setup
    conda uninstall python
    conda install "python=3.10"

Then, build the PolyTop, PolyConf and PolyBuild packages:

.. code-block:: python

    cd polytop

    pip install -e .

    cd ../polyconf

    pip install -e .

    cd ../polybuild

    pip install -e .



Important note: When using PolyTop and PolyBuild (and thus PolyConf to produce matching coordinates), it is always necessary for at least one of the molecules being joined at a polymerization bond to lose at least 2 atoms in its leaving ‘residue’ for sufficient information to exist to infer a dihedral constraint to rotation around the polymerization bond. 

Note that when more than one type of junction exists within the polymer,  it is important that each junction type is given a unique name. In the case where there exist multiple junctions in either molecule sharing the same name, the specific junctions chosen will be randomly distributed among junctions with the same name, allowing for stochastic extension of polymers.  For repeatability it is therefore necessary to use a consistent seed value (in python), and use polytop as a python library rather than from the command line.


Examples:

