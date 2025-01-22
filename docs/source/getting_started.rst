Getting Started
===============

This page details how to get started with PolyConstruct.


**Installing PolyConstruct:**

From your home directory, clone PolyConstruct from Git:

.. code-block:: python

    cd ~
    git clone https://github.com/OMaraLab/polyconstruct.git

Then navigate to polyconstruct:

.. code-block:: python

    cd polyconstruct

To setup an environment for polyconstruct, run: 

.. code-block:: python

    conda create --name polyconstruct-env
    conda activate polyconstruct-env

    pip install -r requirements.txt

    # polytop requires python 3.10
    # you may need to install python=3.10 before you are able to install
    # requirements (as per above), depending on your setup
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


.. note::
    When using PolyTop and PolyBuild (and thus PolyConf to produce matching coordinates),
    it is always necessary for at least one of the molecules being joined at a polymerization
    bond to lose at least 2 atoms in its leaving ‘residue’ for sufficient information to
    exist to infer a dihedral constraint to rotation around the polymerization bond. 

.. note::
    Note that when more than one type of junction exists within a PolyTop polymer,
    it is important that each junction type is given a unique name. In the case where
    there exist multiple junctions in either molecule sharing the same name, the specific
    junctions chosen will be randomly distributed among junctions with the same name,
    allowing for stochastic extension of polymers. For repeatability it is therefore
    necessary to use a consistent seed value (in python), and use PolyTop as a python
    library rather than from the command line. If an exact structure is desired instead,
    simply ensure that each junction type has a unique name that does not allow for any
    discrepancy in exactly which junctions are joined and where.



Worked Examples
==================

**PolyConf**


Simple example - construction of a PEI homopolymer:

.. code-block:: python

    # Import required classes from PolyConf
    from polyconf.Monomer import Monomer
    from polyconf.Polymer import Polymer
    from polyconf.PDB import PDB

    # Initialise Polymer from Monomer of the starting monomer PDB
    polymer=Polymer(Monomer('PEI_start.pdb'))
    imax=127 # define constant, to add an additional 127 monomers

    # Extend the Polymer to the desired length (in this case 128)
    for i in range(0, imax):
        if not i==imax:
            # Extend the Polymer for every step except the last one
            # Extend by one monomer, WITHOUT aligning along this step's linearization vector
            polymer.extend(Monomer('PEI_monomer.pdb'), # extend with this monomer
                n=polymer.maxresid(), # extend existing residue i
                nn=polymer.newresid(), # incoming monomer will have resid i+1
                names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
                joins=[('N1','C1')]) # new connection between N1_i and C1_i+1 
        else:
            # Extend and cap the Polymer by adding the terminating monomer
            polymer.extend(Monomer('PEI_end.pdb'),i,i+1,
                names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'),
                joins=[('N1','C1')])

    # Save the polymer to a file without the dummy atoms, visually check the PDB with another package such as VMD
    Saver = PDB(polymer)
    Saver.cleanup() # center the Polymer in the PBC box
    Saver.save(dummyAtoms='CX NX',fname='polymer_01_vanilla-extend') # save, excluding dummy atoms

    # When you examine the polymer, you can see that the resulting strucure is a tightly coiled helix, rather than linear/
    # This structure is highly ordered, and the turns of the helix are very close.


All of the monomer PDB files used in the above example and the resulting
polymer file are readily available at 'polyconstruct/polyconf_examples/'.



**PolyTop**


Simple example - construction of a linear homopolymer:

.. code-block:: python

    # Import required classes from PolyTop
    from polytop.polytop import Topology, Junction, Monomer, Polymer, Visualize

    # Load in monomer Topology from ITP file
    top = Topology.from_ITP("data_paper_examples/pei.itp")

    # Create a Junction to join 'to' and another to join 'from'.
    # Provide the bonding atom and the leaving atom, in that order, for the
    # Junction - they must have a bond between them.
    to_j = Junction(top.get_atom("C51"), top.get_atom("C62"), name = "to")
    from_j = Junction(top.get_atom("N7"), top.get_atom("C6"), name = "from")

    # Create a Monomer from the Topology and a list of the Junctions
    monomer = Monomer(top, [to_j, from_j])

    # Start the Polymer with one Monomer
    polymer = Polymer(monomer)

    # Extend the Polymer to the desired length (in this case 20)
    for i in range(19):
        polymer.extend(monomer, from_junction_name="from", to_junction_name="to")

    # Save the polymer to a file and visualise the structure with RDKit for an easy visual structure check
    polymer.topology.title = "pei polymer" # renames the ITP header and image
    polymer.save_to_file('data_paper_examples/pei_linear_polymer.json') # text dump
    polymer.topology.to_ITP('data_paper_examples/pei_linear_polymer.itp')
    Visualize.polymer(polymer,infer_bond_order=False).draw2D('data_paper_examples/pei_linear_polymer.png',(400,300))




Complex example - construction of a 4-arm PEG star polymer from single monomeric units:

.. code-block:: python

    # Import required classes from PolyTop
    from polytop.polytop import Topology, Junction, Monomer, Polymer, Visualize

    # Load in monomer topologies from ITP files
    ethanol = Topology.from_ITP("data_paper_examples/extended_ethanol.itp") # main arm monomer
    methane = Topology.from_ITP("data_paper_examples/extended_methane.itp") # terminal monomer
    neopentane = Topology.from_ITP("data_paper_examples/extended_neopentane.itp") # central monomer

    # Create junctions for each monomer with the bonding atom and then the leaving
    # atom specified, in that order, with a unique name. Note how each junction
    # has a unique, descriptive name.
    oxy_j1 = Junction(ethanol.get_atom("O1"), ethanol.get_atom("C1"), name = "oxy1")
    carb_j1 = Junction(ethanol.get_atom("C3"), ethanol.get_atom("O2"), name = "carb1")
    oxy_j2 = Junction(ethanol.get_atom("O1"), ethanol.get_atom("C1"), name = "oxy2")
    carb_j2 = Junction(ethanol.get_atom("C3"), ethanol.get_atom("O2"), name = "carb2")
    oxy_j3 = Junction(ethanol.get_atom("O1"), ethanol.get_atom("C1"), name = "oxy3")
    carb_j3 = Junction(ethanol.get_atom("C3"), ethanol.get_atom("O2"), name = "carb3")
    oxy_j4 = Junction(ethanol.get_atom("O1"), ethanol.get_atom("C1"), name = "oxy4")
    carb_j4 = Junction(ethanol.get_atom("C3"), ethanol.get_atom("O2"), name = "carb4")

    j1 = Junction(neopentane.get_atom("C1"), neopentane.get_atom("O1"), name = "branch1")
    j2 = Junction(neopentane.get_atom("C3"), neopentane.get_atom("O2"), name = "branch2")
    j3 = Junction(neopentane.get_atom("C4"), neopentane.get_atom("O3"), name = "branch3")
    j4 = Junction(neopentane.get_atom("C5"), neopentane.get_atom("O4"), name = "branch4")

    term_j = Junction(methane.get_atom("C1"), methane.get_atom("O1"), name = "term")

    # Create monomers from their topologies and any specified junctions
    e1 = Monomer(ethanol, [oxy_j1, carb_j1])
    e2 = Monomer(ethanol, [oxy_j2, carb_j2])
    e3 = Monomer(ethanol, [oxy_j3, carb_j3])
    e4 = Monomer(ethanol, [oxy_j4, carb_j4])

    central = Monomer(neopentane, [j1, j2, j3, j4])

    terminal = Monomer(methane, [term_j]) # only needs one junction to join to the ends of each arm

    # Start the polymer with the central monomer
    four_polymer = Polymer(central)

    # Attach three ethanols to each of the four junctions (j1-j4) of the central monomer.
    # Note how the extension is done layer by layer, and each of the four arms
    # uses differently named junctions - this ensures that there is no unexpected
    # variation or randomness introduced from multiple degenerately named junctions.
    four_polymer.extend(e1, from_junction_name="branch1", to_junction_name="oxy1")
    four_polymer.extend(e2, from_junction_name="branch2", to_junction_name="oxy2")
    four_polymer.extend(e3, from_junction_name="branch3", to_junction_name="oxy3")
    four_polymer.extend(e4, from_junction_name="branch4", to_junction_name="oxy4")

    four_polymer.extend(e1, from_junction_name="carb1", to_junction_name="oxy1")
    four_polymer.extend(e2, from_junction_name="carb2", to_junction_name="oxy2")
    four_polymer.extend(e3, from_junction_name="carb3", to_junction_name="oxy3")
    four_polymer.extend(e4, from_junction_name="carb4", to_junction_name="oxy4")

    four_polymer.extend(e1, from_junction_name="carb1", to_junction_name="oxy1")
    four_polymer.extend(e2, from_junction_name="carb2", to_junction_name="oxy2")
    four_polymer.extend(e3, from_junction_name="carb3", to_junction_name="oxy3")
    four_polymer.extend(e4, from_junction_name="carb4", to_junction_name="oxy4")

    four_polymer.extend(terminal, from_junction_name="carb1", to_junction_name="term")
    four_polymer.extend(terminal, from_junction_name="carb2", to_junction_name="term")
    four_polymer.extend(terminal, from_junction_name="carb3", to_junction_name="term")
    four_polymer.extend(terminal, from_junction_name="carb4", to_junction_name="term")

    # check polymer charge and give the polymer a descriptive name
    print(f"netcharge = {four_polymer.topology.netcharge}")
    four_polymer.topology.title = "four arm star polymer" # renames the ITP header and image

    # save the polymer to a file and visualise the structure with RDKit for an easy visual structure check
    four_polymer.save_to_file('data_paper_examples/four_arm_star_overlapped_monomers.json') # text dump
    four_polymer.topology.to_ITP('data_paper_examples/four_arm_star_overlapped_monomers.itp')
    Visualize.polymer(four_polymer,infer_bond_order=False).draw2D('data_paper_examples/four_arm_star_overlapped_monomers.png',(400,300))

All of the monomer ITP files used in the above two examples, and the resulting
polymer files, are also readily available at 'polyconstruct/data_paper_examples/'.


**Automated Builder**

An automated builder is available from the command line, to simultaneously
construct united-atom, linear polymer structures and parameters. 

This builder streamlines the process of making matching polymer structure and
parameter files, as it leverages PolyConf and PolyTop to build the polymer
simultaneously. Additionally, no coding is required to use this module, as it
is fully available from the command line.

This automated builder is recently developed and still undergoing development
to ensure it is robust and suitable for all cases. Please manually check the
topology and structure files produced by this module for correctness (e.g. with
GROMACS 'gmx grompp' and visual inspection with VMD).

The builder is not working for your case? Please open an issue on the
`PolyConstruct GitHub repository <https://github.com/OMaraLab/polyconstruct>`_
with the command used, a copy of your monomer CSV file and a description of the
issue observed.

For more information and how to use it, see
:ref:`Automated Linear, United-Atom Polymer Builder` documentation.


--------------------------------------------------------------------------------------

Find the above and additional worked examples as executable Python scripts or Jupyter
Notebooks on the `PolyConstruct GitHub repository <https://github.com/OMaraLab/polyconstruct>`_.
Examples for PolyTop are available at 'polyconstruct/paper_worked_examples.ipynb'
and for PolyConf at 'polyconstruct/polyconf_examples/'; while instructions to
use the two PolyBuild scripts are included under the :ref:`PolyBuild` documentation.