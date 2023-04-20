# Polytop Library Demonstration

This Jupyter notebook will demonstrate how to use the Polytop library to load sample ITP files, create monomer instances, and display a 3D representation of molecules.

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
    



We can render the topology using 3Dmol (py3DMol) to see what it looks like.


```python
# render display presentation of arginine
import py3Dmol
from polytop.visualize import Visualize
viewer = py3Dmol.view(width=800, height=400)
Visualize(arg).create_py3Dmol_view(viewer)
viewer.show()
```


<div id="3dmolviewer_16820058569475915"  style="position: relative; width: 800px; height: 400px">
        <p id="3dmolwarning_16820058569475915" style="background-color:#ffcccc;color:black">You appear to be running in JupyterLab (or JavaScript failed to load for some other reason).  You need to install the 3dmol extension: <br>
        <tt>jupyter labextension install jupyterlab_3dmol</tt></p>
        </div>
<script>

var loadScriptAsync = function(uri){
  return new Promise((resolve, reject) => {
    //this is to ignore the existence of requirejs amd
    var savedexports, savedmodule;
    if (typeof exports !== 'undefined') savedexports = exports;
    else exports = {}
    if (typeof module !== 'undefined') savedmodule = module;
    else module = {}

    var tag = document.createElement('script');
    tag.src = uri;
    tag.async = true;
    tag.onload = () => {
        exports = savedexports;
        module = savedmodule;
        resolve();
    };
  var firstScriptTag = document.getElementsByTagName('script')[0];
  firstScriptTag.parentNode.insertBefore(tag, firstScriptTag);
});
};

if(typeof $3Dmolpromise === 'undefined') {
$3Dmolpromise = null;
  $3Dmolpromise = loadScriptAsync('https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.1/3Dmol-min.js');
}

var viewer_16820058569475915 = null;
var warn = document.getElementById("3dmolwarning_16820058569475915");
if(warn) {
    warn.parentNode.removeChild(warn);
}
$3Dmolpromise.then(function() {
viewer_16820058569475915 = $3Dmol.createViewer(document.getElementById("3dmolviewer_16820058569475915"),{backgroundColor:"white"});
viewer_16820058569475915.zoomTo();
	viewer_16820058569475915.addModel("\n     RDKit          3D\n\n 26 25  0  0  0  0  0  0  0  0999 V2000\n   -4.9951    1.5338   -0.6108 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.2088    1.7752   -0.0201 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.7249    2.5864   -0.3951 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.3359    0.7532    0.1060 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.9971   -0.3852    0.4094 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.4059   -1.1235    0.7680 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7767   -0.2138    1.0344 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0865    0.9576   -0.0677 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1740   -0.1737    0.0520 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1766   -0.5585    1.0788 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4668   -0.9819   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.2313    0.3102   -0.3067 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5103    1.1334    0.3618 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.2291    0.7174   -1.3256 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2538   -0.8255   -0.2055 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9192   -1.6527   -0.8452 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2689   -1.2056    0.8244 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6622   -0.3908   -0.6317 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6439    0.0593   -1.6301 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.5625   -1.5671   -0.6878 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7368   -1.8712    0.2733 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0918   -2.3381   -1.1566 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2357    0.6437    0.3439 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2668    1.8537    0.1864 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.6983    0.0947    1.4923 O   0  0  0  0  0  0  0  0  0  0  0  0\n    4.0376    0.8691    1.9878 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  4  1  0\n  2  3  1  0\n  4  8  2  3\n  4  5  1  0\n  5  7  1  0\n  5  6  1  0\n  8  9  1  0\n  9 12  1  0\n  9 11  1  0\n  9 10  1  0\n 12 14  1  0\n 12 13  1  0\n 12 15  1  0\n 15 17  1  0\n 15 18  1  0\n 15 16  1  0\n 18 23  1  0\n 18 19  1  0\n 18 20  1  0\n 20 22  1  0\n 20 21  1  0\n 23 25  1  0\n 23 24  2  0\n 25 26  1  0\nM  END\n","mol");
	viewer_16820058569475915.setStyle({"stick": {}});
	viewer_16820058569475915.zoomTo();
viewer_16820058569475915.render();
});
</script>


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
glu_monomer = Monomer(glu,start,end).save('tests/samples/glutamine_monomer.json')

Visualize(arg).infer_bond_orders().create_2D_image('tests/samples/arginine.png',(200,100))
Visualize(arg_monomer.LHS).infer_bond_orders().create_2D_image('tests/samples/arginine_LHS.png',(200,100))
Visualize(arg_monomer.link).infer_bond_orders().create_2D_image('tests/samples/arginine_link.png',(200,100))
Visualize(arg_monomer.RHS).infer_bond_orders().create_2D_image('tests/samples/arginine_RHS.png',(200,100))
from IPython.display import Image
from IPython.core.display import HTML
html = ''
html += f'<figure><img src="tests/samples/arginine.png" style="margin:0 10px" width="300"><figcaption>Arginine molecule</figcaption></figure>'

arginines = ['arginine_LHS.png','arginine_link.png', 'arginine_RHS.png']
html += '<div style="display:flex">'
for image in arginines:
    html += f'<figure><img src="tests/samples/{image}" style="margin:0 10px" width="300"><figcaption>{image}</figcaption></figure>'
html += '</div>'
display(HTML(html))

```


<figure><img src="tests/samples/arginine.png" style="margin:0 10px" width="300"><figcaption>Arginine molecule</figcaption></figure><div style="display:flex"><figure><img src="tests/samples/arginine_LHS.png" style="margin:0 10px" width="300"><figcaption>arginine_LHS.png</figcaption></figure><figure><img src="tests/samples/arginine_link.png" style="margin:0 10px" width="300"><figcaption>arginine_link.png</figcaption></figure><figure><img src="tests/samples/arginine_RHS.png" style="margin:0 10px" width="300"><figcaption>arginine_RHS.png</figcaption></figure></div>


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
