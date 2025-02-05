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
conda create --name polyconstruct-env
conda activate polyconstruct-env

pip install -r requirements.txt
```

If you recieve an error in the last step, it may be from the Python version on
your computer setup, as PolyTop has a strict requirement for version 3.10. Run
the below commands then repeat `pip install -r requirements.txt` :

```
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
