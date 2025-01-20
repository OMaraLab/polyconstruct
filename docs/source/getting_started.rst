Getting Started
===============

This page details how to get started with PolyConstruct.

Important note: When using PolyTop and PolyBuild (and thus PolyConf to produce matching coordinates), it is always necessary for at least one of the molecules being joined at a polymerization bond to lose at least 2 atoms in its leaving ‘residue’ for sufficient information to exist to infer a dihedral constraint to rotation around the polymerization bond. 

Note that when more than one type of junction exists within the polymer,  it is important that each junction type is given a unique name. In the case where there exist multiple junctions in either molecule sharing the same name, the specific junctions chosen will be randomly distributed among junctions with the same name, allowing for stochastic extension of polymers.  For repeatability it is therefore necessary to use a consistent seed value (in python), and use polytop as a python library rather than from the command line.  