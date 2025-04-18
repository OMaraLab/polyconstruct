# PolyTop example scripts

This folder contains a handful of example scripts demonstrating the use of
*PolyTop*. All the files you need to run these examples are contained within
this folder.

For detailed documentation of all *PolyTop* methods, please refer to the
[PolyConstruct ReadTheDocs](https://polyconstruct.readthedocs.io/en/latest/index.html).

These tutorials are not intended as a primer on polymer design or
stereochemistry. It is assumed that the reader is familiar with polymer
chemistry and common types of copolymers.

## Tutorial 01: Building a linear polyethyleneimine polymer topology using PolyTop

The first tutorial is an example of how to make a linear 20 unit
polyethyleneimine chain using *PolyTop*. It is contained in the file
`linear_PEI.py`.

This tutorial showcases several approaches, including:
* loading in itp files to create Topology objects
* defining appropriate Junctions for polymer extension
* creating Monomer objects from a Topology and list of Junctions
* building a polymer chain with `extend()`
* saving the polymer to an itp file

This tutorial is all contained in a single script, with detailed comments
discussing the strategies used and the reasoning behind them or another method.
You are encouraged to view the topologies generated by this script, and others
you develop, to understand how the monomers are connected.

## Tutorial 02: Building an ethlyamine dendrimer polymer topology using PolyTop

The second tutorial shows strategies for making a more complex polymer with
levels of branching building outwards from a central monomer. It is contained
in the file `dendrimer_ethylamine.py`.

This tutorial showcases several approaches, including:
* breaking down a polymer structure into the monomers required to build it
* loading in itp files to create Topology objects
* defining appropriate, distinct Junctions for precise polymer extension
* creating Monomer objects from a Topology and list of Junctions
* building a polymer chain with `extend()`
* saving the polymer to an itp file

This tutorial is all contained in a single script, with detailed comments
discussing the strategies used and the reasoning behind them or another method.
You are encouraged to view the topologies generated by this script, and others
you develop, to understand how the monomers are connected.

## Tutorial 03: Building a 4-arm PEG star polymer topology using PolyTop

The third tutorial shows strategies for making another complex polymer with
long linear arms branching out from a central monomer. It is contained
in the file `star_PEG.py`.

This tutorial showcases several approaches, including:
* breaking down a polymer structure into the monomers required to build it
* discussing different approaches to create the same polymer structure
* loading in itp files to create Topology objects
* defining appropriate, distinct Junctions for precise polymer extension
* creating Monomer objects from a Topology and list of Junctions
* building a polymer chain with `extend()`
* saving the polymer to an itp file

This tutorial is all contained in a single script, with detailed comments
discussing the strategies used and the reasoning behind them or another method.
You are encouraged to view the topologies generated by this script, and others
you develop, to understand how the monomers are connected.
