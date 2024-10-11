Setup for Polytop

From your home directory, install polytop from Git:

```
git clone https://github.com/OMaraLab/polytop.git
```

Then naviate to polytop:

```
cd polytop
```

To setup polytop, run: 

```
conda create --name polytop-env
conda activate polytop-env

pip install -r requirements.txt

# polytop requires python 3.10, and setup.py will not run without it
conda uninstall python
conda install "python=3.10" 

cd polytop

pip install -e .

cd ../polyconf

pip install -e .
```
