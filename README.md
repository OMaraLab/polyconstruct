# PolyConstruct: a python library for polymer generation

PolyConstruct contains three python tools for generating polymer coordinate and topolgy files for molecular dynamics simulations.  It was published in the paper *PolyConstruct: adapting  biomolecular simulation pipelines for polymers with PolyBuild, PolyConf and PolyTop*.

## PolyConstruct guides

Detailed documentation is available at the [PolyConstruct ReadTheDocs](https://polyconstruct.readthedocs.io/en/latest/index.html).  This includes installation instructions, tutorials and worked examples, and api documentation for all polyconstruct methods.

There are also detailed tutorials for PolyConf in the folder [polyconf_examples](https://github.com/OMaraLab/polyconstruct/tree/main/polyconf_examples), and detailed tutorials for PolyTop in the folder [polytop_examples](https://github.com/OMaraLab/polyconstruct/tree/main/polytop_examples)

## Quick and dirty installation

From your home directory, clone PolyConstruct from Git:

```bash
cd ~
git clone https://github.com/OMaraLab/polyconstruct.git
```

Then navigate to polyconstruct:

```bash
cd polyconstruct
```

To setup polyconstruct, run:

```bash
conda create --name polyconstruct python=3.10
conda activate polyconstruct

pip install -r requirements.txt
```

Then, build the PolyTop, PolyConf and PolyBuild packages:

```bash
cd polytop

pip install -e .

cd ../polyconf

pip install -e .

cd ../polybuild

pip install -e .
```
