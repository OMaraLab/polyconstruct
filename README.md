# PolyConstruct:  adapting  biomolecular simulation pipelines for polymers with PolyBuild, PolyConf and PolyTop 

## Documentation

The full documentation is available at the [PolyConstruct ReadTheDocs](https://polyconstruct.readthedocs.io/en/latest/index.html).  This includes installation instructions, documentation for all polyconstruct methods, and several worked examples and tutorials.

There are also example scripts demonstrating the use of PolyConf in the folder [polyconf_examples](https://github.com/OMaraLab/polyconstruct/tree/main/polyconf_examples), and example scripts demonstrating the use of PolyTop in the folder [polytop_examples](https://github.com/OMaraLab/polyconstruct/tree/main/polytop_examples)

## How to install PolyConstruct

From your home directory, install PolyConstruct from Git:

```
cd ~
git clone https://github.com/OMaraLab/polyconstruct.git
```

Then navigate to polyconstruct:

```
cd polyconstruct
```

To setup polyconstruct, run: 

```
conda create --name polyconstruct python=3.10
conda activate polyconstruct

pip install -r requirements.txt
```

Then, build the PolyTop, PolyConf and PolyBuild packages:

```
cd polytop

pip install -e .

cd ../polyconf

pip install -e .

cd ../polybuild

pip install -e .
```

