# Polytop Library Demonstration

This Jupyter notebook will demonstrate how to use the Polytop library to load sample ITP files, create monomer instances, and display a 2D representation of molecules.

## Importing the Library

We will import the `polytop` library and confirm the version.


```python
# Load sample ITPs
import polytop
from polytop.monomer import Monomer 
from polytop.topology import Topology
from polytop.visualize import Visualize

print(f"Polytop Version {polytop.__version__}")
```

    Polytop Version 0.1 alpha
    

## Loading Sample ITP Files

First, we will load the ITP files for some molecules to use as monomers.  We'll use some amino acids to make a polypeptide.  We'll start by importing `Arginine` using the `Topology.from_ITP()` function. This function will return a `Topology` object that contains all of the information from the ITP file.  Note for display purposes we'll add a double bond manually between the carbon and the second oxygen in the terminal carboxyl group.


```python
arg = Topology.from_ITP('tests/samples/arginine.itp')
```

We can visualize the topology as a 2D structure using the Visualize class.


```python
Visualize(arg).create_2D_image('README_files/arginine.png',(400,200))
from IPython.display import Image
Image(filename='README_files/arginine.png') 
```




    
![png](README_files/README_5_0.png)
    



NB: ITP topologies do not contain information about bond ordering.  We can manually fix up any bond.


```python
arg.get_bond(23,24).order = 2 # ITP files do not include bond orders
arg.get_bond(4,8).order = 2
Visualize(arg).create_2D_image('README_files/arginine.png',(400,200))
from IPython.display import Image
Image(filename='README_files/arginine.png')
```




    
![png](README_files/README_7_0.png)
    



We can also infer bond orders from expected atom valencies and available bonds.  This is done by the Visualize.infer_bond_orders() function.


```python
arg = Topology.from_ITP('tests/samples/arginine.itp')
Visualize(arg).infer_bond_orders().create_2D_image('README_files/arginine.png',(400,200))
from IPython.display import Image
Image(filename='README_files/arginine.png')
```




    
![png](README_files/README_9_0.png)
    



## Loading a glutamine topology

And we'll do the same with the second amino acid, glutamine


```python
glu = Topology.from_ITP('tests/samples/glutamine.itp')
Visualize(glu).infer_bond_orders().create_2D_image('README_files/glutamine.png',(400,200))

from IPython.display import Image
Image(filename='README_files/glutamine.png')
```




    
![png](README_files/README_11_0.png)
    



# Convert a molecule topology to a monomer 

A monomer is an element that can participate in the formation of a polymer.  A monomer is defined by a topology and a set of polymerization junctions (defined by start and end bonds).  The monomer has 3 potential topologies, a start  unit, a link unit, and an end unit.  Which variant will be used will be determined by the position in the polymer the unit occupies.  The start unit is used for the first monomer in the polymer, the end unit is used for the last monomer in the polymer, and the link unit is used for all other monomers. 


```python
arg = Topology.from_ITP('tests/samples/arginine.itp')
start = arg.get_bond(21,20)
end = arg.get_bond(23,25)
arg_monomer = Monomer(arg,start,end)
arg_monomer.save('tests/samples/arginine_monomer.json')

glu=Topology.from_ITP('tests/samples/glutamine.itp')
start = glu.get_bond(7,8)
end = glu.get_bond(2,3)
glu_monomer = Monomer(glu,start,end)
glu_monomer.save('tests/samples/glutamine_monomer.json')

Visualize(arg).infer_bond_orders().create_2D_image('tests/samples/arginine.png',(400,200))
Visualize(arg_monomer.LHS).infer_bond_orders().create_2D_image('tests/samples/arginine_LHS.png',(400,200))
Visualize(arg_monomer.link).infer_bond_orders().create_2D_image('tests/samples/arginine_link.png',(400,200))
Visualize(arg_monomer.RHS).infer_bond_orders().create_2D_image('tests/samples/arginine_RHS.png',(400,200))

Visualize(glu).infer_bond_orders().create_2D_image('tests/samples/glutamine.png',(400,200))
Visualize(glu_monomer.LHS).infer_bond_orders().create_2D_image('tests/samples/glutamine_LHS.png',(400,200))
Visualize(glu_monomer.link).infer_bond_orders().create_2D_image('tests/samples/glutamine_link.png',(400,200))
Visualize(glu_monomer.RHS).infer_bond_orders().create_2D_image('tests/samples/glutamine_RHS.png',(400,200))

from IPython.display import Image
from IPython.core.display import HTML
html = '<h1>Arginine</h1>'
html += f'<figure><img src="tests/samples/arginine.png" style="margin:0 10px" width="300"><figcaption>Arginine molecule</figcaption></figure>'

arginines = ['arginine_LHS.png','arginine_link.png', 'arginine_RHS.png']
html += '<div style="display:flex">'
for image in arginines:
    html += f'<div style="display:inline-flex"><figure><img src="tests/samples/{image}" style="margin:0 10px" width="300"><figcaption>{image}</figcaption></figure></div>'
html += '</div>'

html += '<h1>Glutamine</h1>'
html += f'<figure><img src="tests/samples/glutamine.png" style="margin:0 10px" width="300"><figcaption>Glutamine molecule</figcaption></figure>'

glutamines = ['glutamine_LHS.png','glutamine_link.png', 'glutamine_RHS.png']
html += '<div style="display:flex">'
for image in glutamines:
    html += f'<div style="display:inline-flex"><figure><img src="tests/samples/{image}" style="margin:0 10px" width="300"><figcaption>{image}</figcaption></figure></div>'
html += '</div>'

display(HTML(html))

```


<h1>Arginine</h1><figure><img src="tests/samples/arginine.png" style="margin:0 10px" width="300"><figcaption>Arginine molecule</figcaption></figure><div style="display:flex"><div style="display:inline-flex"><figure><img src="tests/samples/arginine_LHS.png" style="margin:0 10px" width="300"><figcaption>arginine_LHS.png</figcaption></figure></div><div style="display:inline-flex"><figure><img src="tests/samples/arginine_link.png" style="margin:0 10px" width="300"><figcaption>arginine_link.png</figcaption></figure></div><div style="display:inline-flex"><figure><img src="tests/samples/arginine_RHS.png" style="margin:0 10px" width="300"><figcaption>arginine_RHS.png</figcaption></figure></div></div><h1>Glutamine</h1><figure><img src="tests/samples/glutamine.png" style="margin:0 10px" width="300"><figcaption>Glutamine molecule</figcaption></figure><div style="display:flex"><div style="display:inline-flex"><figure><img src="tests/samples/glutamine_LHS.png" style="margin:0 10px" width="300"><figcaption>glutamine_LHS.png</figcaption></figure></div><div style="display:inline-flex"><figure><img src="tests/samples/glutamine_link.png" style="margin:0 10px" width="300"><figcaption>glutamine_link.png</figcaption></figure></div><div style="display:inline-flex"><figure><img src="tests/samples/glutamine_RHS.png" style="margin:0 10px" width="300"><figcaption>glutamine_RHS.png</figcaption></figure></div></div>


# Convert monomers + distribution to a Polymer

Once you have multiple molecules configured as monomers you can then polymerize them.  The polymer will be extended with a specific number of monomer units (here 12), with a specific random uniform distribution (20% arginine / 80% glutamine), and optionally using a random seed (42) and optionally a specific start monomer and end monomer.  The polymer can be saved to a .json file 


```python

from polytop.polymer import Polymer


polymer = Polymer([arg_monomer,glu_monomer], [20,80], num_monomers= 12, seed= 42, start_monomer= arg_monomer)
polymer.save_to_file('tests/samples/polymer.json')
polymer_topology = polymer.get_topology()

Visualize(polymer_topology).infer_bond_orders().create_2D_image('tests/samples/polymer.png',(400,200))
from IPython.display import Image
Image(filename='README_files/polymer.png')
```


    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    c:\Users\Richard\OneDrive - Australian National University\Polytop\polytop\README.ipynb Cell 18 in <cell line: 5>()
          <a href='vscode-notebook-cell:/c%3A/Users/Richard/OneDrive%20-%20Australian%20National%20University/Polytop/polytop/README.ipynb#X23sZmlsZQ%3D%3D?line=0'>1</a> from polytop.polymer import Polymer
          <a href='vscode-notebook-cell:/c%3A/Users/Richard/OneDrive%20-%20Australian%20National%20University/Polytop/polytop/README.ipynb#X23sZmlsZQ%3D%3D?line=3'>4</a> polymer = Polymer([arg_monomer,glu_monomer], [20,80], num_monomers= 12, seed= 42, start_monomer= arg_monomer)
    ----> <a href='vscode-notebook-cell:/c%3A/Users/Richard/OneDrive%20-%20Australian%20National%20University/Polytop/polytop/README.ipynb#X23sZmlsZQ%3D%3D?line=4'>5</a> polymer.save_to_file('tests/samples/polymer.json')
          <a href='vscode-notebook-cell:/c%3A/Users/Richard/OneDrive%20-%20Australian%20National%20University/Polytop/polytop/README.ipynb#X23sZmlsZQ%3D%3D?line=5'>6</a> polymer_topology = polymer.get_topology()
          <a href='vscode-notebook-cell:/c%3A/Users/Richard/OneDrive%20-%20Australian%20National%20University/Polytop/polytop/README.ipynb#X23sZmlsZQ%3D%3D?line=7'>8</a> Visualize(polymer_topology).infer_bond_orders().create_2D_image('tests/samples/polymer.png',(400,200))
    

    File c:\Users\Richard\OneDrive - Australian National University\Polytop\polytop\polytop\polymer.py:58, in Polymer.save_to_file(self, filename)
         56 def save_to_file(self, filename: str) -> None:
         57     with open(filename, "w") as f:
    ---> 58         json.dump(self.to_dict(), f)
    

    File c:\Users\Richard\OneDrive - Australian National University\Polytop\polytop\polytop\polymer.py:76, in Polymer.to_dict(self)
         65 def to_dict(
         66     self,
         67 ) -> Dict[
       (...)
         73     ],
         74 ]:
         75     return {
    ---> 76         "monomers": [monomer.to_dict() for monomer in self.monomers],
         77         "distribution": self.distribution,
         78         "num_monomers": self.num_monomers,
         79         "seed": self.seed,
         80         "start_monomer": self.start_monomer.to_dict()
         81         if self.start_monomer
         82         else None,
         83         "end_monomer": self.end_monomer.to_dict() if self.end_monomer else None,
         84     }
    

    File c:\Users\Richard\OneDrive - Australian National University\Polytop\polytop\polytop\polymer.py:76, in <listcomp>(.0)
         65 def to_dict(
         66     self,
         67 ) -> Dict[
       (...)
         73     ],
         74 ]:
         75     return {
    ---> 76         "monomers": [monomer.to_dict() for monomer in self.monomers],
         77         "distribution": self.distribution,
         78         "num_monomers": self.num_monomers,
         79         "seed": self.seed,
         80         "start_monomer": self.start_monomer.to_dict()
         81         if self.start_monomer
         82         else None,
         83         "end_monomer": self.end_monomer.to_dict() if self.end_monomer else None,
         84     }
    

    AttributeError: 'NoneType' object has no attribute 'to_dict'


To convert this notebook to markdown, run the following command: 
    ```
    jupyter nbconvert README.ipynb --to markdown
    ```


```python

```
