# *PolyConstruct*: a python library for polymer generation

*PolyConstruct* contains three python tools for generating polymer coordinate and topolgy files for molecular dynamics simulations.  

* *PolyConf* is a tool for generating ensembles of polymer conformations by combining monomer coordinate files.

* *PolyBuild* is a tool for generating polymer topology files for simulaton from polymer coordinate files, leveraging the functionality of the gromacs tool pdb2gmx.

* *PolyTop* is a tool for generating polymer topology files from monomer topology files.

*PolyConstruct* was published in the paper *PolyConstruct: adapting biomolecular simulation pipelines for polymers with PolyBuild, PolyConf and PolyTop*.

## Getting started

Detailed documentation is available at the [PolyConstruct ReadTheDocs](https://polyconstruct.readthedocs.io/en/latest/index.html).  This includes installation instructions, tutorials and worked examples, and api documentation for all *PolyConstruct* methods.

There is a series of detailed tutorials for *PolyConf* in the folder [polyconf_examples](https://github.com/OMaraLab/polyconstruct/tree/main/polyconf_examples), and for *PolyTop* in the folder [polytop_examples](https://github.com/OMaraLab/polyconstruct/tree/main/polytop_examples).  A set of example input and output files for *PolyBuild* are presented in the folder [polybuild_examples/RTP_entries](https://github.com/OMaraLab/polyconstruct/tree/main/polybuild_examples/RTP_entries)

## Quick and dirty installation

From your home directory, clone *PolyConstruct* from this repository:

```bash
cd ~
git clone https://github.com/OMaraLab/polyconstruct.git
```

Then navigate to `~/polyconstruct`

```bash
cd polyconstruct
```

Create a python environment and setup *PolyConstruct*:

```bash
conda create --name polyconstruct python=3.10
conda activate polyconstruct

pip install -r requirements.txt
```

Then, build the *PolyTop*, *PolyConf* and *PolyBuild* packages:

```bash
cd polytop

pip install -e .

cd ../polyconf

pip install -e .

cd ../polybuild

pip install -e .
```
