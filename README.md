# Setup for PolyConstruct

From your home directory, install PolyConstruct from Git:

```
git clone https://github.com/OMaraLab/polyconstruct.git
```

Then navigate to polyconstruct:

```
cd polyconstruct
```

To setup polyconstruct, run: 

```
conda create --name polyconstruct-env
conda activate polyconstruct-env

pip install -r requirements.txt

# polytop requires python 3.10, and setup.py will not run without it
conda uninstall python
conda install "python=3.10" 

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
