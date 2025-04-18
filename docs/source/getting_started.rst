Getting Started with PolyConstruct
===============

This page details how to get started with *PolyConstruct*.

The design, testing, and validation of *PolyConstruct* is detailed in the manuscript *"PolyConstruct: 
adapting  biomolecular simulation pipelines for polymers with PolyBuild, PolyConf and PolyTop"*, by 
*Rangika Munaweera, Ada Quinn, Luna Morrow, Richard A Morris, and Megan L O’Mara*


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

    conda create --name polyconstruct-env python=3.10
    conda activate polyconstruct-env

    pip install -r requirements.txt

Then, build the PolyTop, PolyConf and PolyBuild packages:

.. code-block:: python

    cd polytop

    pip install -e .

    cd ../polyconf

    pip install -e .

    cd ../polybuild

    pip install -e .



Input files
===============

**PolyConf**

*PolyConf* assembles polymer pdb files by tiling individual monomer pdb files.  Monomer PDB files should 
include all monomer atoms that will be present in the final polymer, dummy atoms corresponding to 
connectivity with adjacent monomers, and `CONECT` records that describe the connectivity between
monomer atoms.

*PolyConf* does not require any preprocessing of monomer .pdb files.

**PolyBuild**

*PolyBuild* requires monomer itp files, which are then converted into gromacs residue topology database 
entries, and a coordinate file containing a complete structure of the intended polymer.  This structure can 
be created with PolyConf, manually, or by using other structure creation packages.

PolyBuild leverages GROMACS tools to automatically parameterise polymers. and requires 
.itp files to be converted to RTP database entries using the provided `ITP2RTP` program.  
For more infomation, see the :ref:`PolyBuild` documentation.

**PolyTop**

*PolyTop* requires monomer itp files describing all monomers that will be contained in the system.

*PolyTop* does not require any preprocessing of monomer .itp files.

Monomer design
===============

For both *PolyBuild* and *PolyTop*, monomer parameters should include all monomer atoms that will be present 
in the simulation model of the final polymer, and dummy atoms correesponding to connectivity with adjacent monomers.  
The method for preparing monomer .itp depends on the choice of force field, and a number of suitable automated tools 
exist for small molecule parameterization such as the `Automated Topology Builder <https://atb.uq.edu.au>`_, 
`antechamber <https://ambermd.org/antechamber/ac.html>`_, and `LigParGen <https://zarbi.chem.yale.edu/ligpargen/>`_. 

Monomer coordinates and parameters used as inputs should be designed to reflect the state of the monomer 
in the mature polymer chain, rather than the isolated monomer molecule prior to polymerization.  For example, 
when parameterising a peptide polymer synthesised through amide condensation reactions, monomers should be 
parameterized with their connectivity described as the amide bonds that will occur in the mature polymer, 
and not the ammonium and carboxylate groups found in the precursor amino acids.

Monomer parameters should be designed in a manner consistent with the desired force field.  We recommend you 
do not combine monomer parameters from different force field families in a single polymer. 


For *PolyConf* These coordinate files should represent a sensible monomer geometry as could be found in the final 
polymer.  These might be created using tools like ChimeraX, with theoretical ideal bond lengths and 
angles, by geometry optimization, or generated by automated parameterization tools.  Additionally, 
it is often convenient if the monomer is in a conformation where atoms corresponding to adjacent monomers 
are as far as possible from other atoms, and pointed away from other attachment points.

Worked Examples
==================

**PolyConf**

PolyConf creates polymer coordinate files through the tiling and manipulation of monomer pdb files.

There are several detailed examples of the use of PolyConf to create ensembles of starting 
conformations for a series of increasingly complex polymer architectures.  These are contained in the PolyConf repository in the folder 'polyconstruct/polyconf_examples/' contained in the 

Here is one simple example, showing the construction of a linear polyethylenimine 128-mer.

.. code-block:: python

    # start by loading PolyConf

    from polyconf.Monomer import Monomer
    from polyconf.Polymer import Polymer
    from polyconf.PDB import PDB

    PEI=Polymer(Monomer('PEI_start.pdb')) # initialize polymer with the starting monomer

    for _ in range (0,126): 	# extend our polymer by adding new monomers 

        n= PEI.maxresid() 	# we will extend the monomer with the highest resid
        new_n= PEI.newresid() 	# generate a new resid for the new monomer

        PEI.extend( 					# extend with one monomer
                Monomer('PEI_monomer.pdb'), 	# create new monomer from pdb file
                n=n,					# resid to extend from
                nn=new_n,				# resid for new monomer
                names=dict(P='C1', Q='CX', R='NX', S='N1'), 	# defines the mapping
                joins=[('N1','C1')], 				# defines the new connectivity
                ) 
        
        # after adding each monomer, we will flag dummy atoms
        PEI.renamer(n,'CX') # convert atom CXn to a dummy atom 
        PEI.renamer(new_n,'NX') # convert atom NXnew_n to a dummy atom

    # now we add the final monomer
    # the process is the same, but we are using a different pdb file 

    n= PEI.maxresid()
    new_n= PEI.newresid()
     
    PEI.extend(
                Monomer('PEI_end.pdb'),
                n=n, 
                nn=new_n,
                names=dict(P='C1',Q='CX',R='NX',S='N1'),
                joins=[('N1','C1')],
                )

    # remove dummy atoms as before

    PEI.renamer(n,'CX') 
    PEI.renamer(new_n,'NX') 

    # generate an ensemble of conformations
    # begin by generating lists of dihedrals

    NC_dh = PEI.gen_pairlist(J='N1',K='C1',first_resid=1,last_resid=127,mult=3,same_res=False) 
    alkane_dh = PEI.gen_pairlist(J='C1',K='C2',first_resid=1,last_resid=128,mult=3)

    # Then, we will generate five starting conformations.  

    # For each conformation:
    #   start by making a copy of our initial conformation
    #   randomize the conformation by shuffling the N1n-C1n+1¬¬ torsions
    #   solve the conformation by rotating the C1n-C2n torsions

    for conf in range(1,6): 

        newconf= PEI.copy() # make a duplicate of the original polymer

        newconf.shuffler(NC_dh)
        newconf.dihedral_solver(alkane_dh,cutoff=0.9)

        Saver = PDB(newconf)
        Saver.cleanup() # place the polymer in the center of the simulation box
        Saver.save(fname=f'PEI_linear_conformation_{conf}')

    # end of example script


**PolyBuild**

Example input itp files and the resulting rtp database entries are presented in the folder `polybuild_examples/RTP_entries` 

RTP entries 

**PolyTop**



Simple example - construction of a linear homopolymer:



.. code-block:: python

    # Import required classes from PolyTop
    from polytop.Junction import Junction
    from polytop.Monomer import Monomer
    from polytop.Visualize import Visualize
    from polytop.Polymer import Polymer
    from polytop.Topology import Topology

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


.. note::
    Note that when more than one type of junction exists within a PolyTop polymer,
    it is important that each junction type is given a unique name. 


Complex example - construction of a 4-arm PEG star polymer from single monomeric units:

.. code-block:: python

    # Import required classes from PolyTop
    from polytop.Junction import Junction
    from polytop.Monomer import Monomer
    from polytop.Visualize import Visualize
    from polytop.Polymer import Polymer
    from polytop.Topology import Topology

    # Load in monomer topologies from ITP files
    ethanol = Topology.from_ITP("polytop_examples/data/extended_ethanol.itp") # main arm monomer
    methane = Topology.from_ITP("polytop_examples/data/extended_methane.itp") # terminal monomer
    neopentane = Topology.from_ITP("polytop_examples/data/extended_neopentane.itp") # central monomer

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
    four_polymer.save_to_file('polytop_examples/data/four_arm_star.json') # text dump
    four_polymer.topology.to_ITP('polytop_examples/data/four_arm_star.itp')
    Visualize.polymer(four_polymer,infer_bond_order=False).draw2D('polytop_examples/data/four_arm_star.png',(800,600))

--------------------------------------------------------------------------------------

Find the above and additional worked examples as executable Python scripts on the `PolyConstruct GitHub repository <https://github.com/OMaraLab/polyconstruct>`_.

Examples for *PolyTop* are available at `polyconstruct/polytop_examples/`

Examples for *PolyConf* at `polyconstruct/polyconf_examples/`

Instructions to use *PolyBuild* are included under the :ref:`PolyBuild` documentation.
