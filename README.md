# PolyConstruct Documentation

[PolyConstruct ReadTheDocs](https://polyconstruct.readthedocs.io/en/latest/).


# Setup for PolyConstruct

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
conda create --name polyconstruct-env python=3.10
conda activate polyconstruct-env

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
