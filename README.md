# Polytop Library Demonstration

This Jupyter notebook will demonstrate how to use the Polytop library to load sample ITP files, create monomers from molecules, create polymers from monomers and a distribution, and to display 2D and 3D representation of the toplogies.

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


<div id="3dmolviewer_16816845087131624"  style="position: relative; width: 800px; height: 400px">
        <p id="3dmolwarning_16816845087131624" style="background-color:#ffcccc;color:black">You appear to be running in JupyterLab (or JavaScript failed to load for some other reason).  You need to install the 3dmol extension: <br>
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

var viewer_16816845087131624 = null;
var warn = document.getElementById("3dmolwarning_16816845087131624");
if(warn) {
    warn.parentNode.removeChild(warn);
}
$3Dmolpromise.then(function() {
viewer_16816845087131624 = $3Dmol.createViewer(document.getElementById("3dmolviewer_16816845087131624"),{backgroundColor:"white"});
viewer_16816845087131624.zoomTo();
	viewer_16816845087131624.addModel("\n     RDKit          3D\n\n 26 25  0  0  0  0  0  0  0  0999 V2000\n   -4.4272    2.1021   -0.0189 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.5795    1.8672    0.4829 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8884    2.6032    0.3676 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.0563    0.6917    0.0745 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.0211   -0.2526    0.0277 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.6719   -1.2007   -0.0223 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.6859   -0.1533    0.7866 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8120    0.6113   -0.2054 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2844   -0.6689   -0.6641 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3722   -1.4299    0.1203 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8342   -1.0136   -1.5481 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.1861   -0.4915   -1.0472 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.2779    0.3088   -1.7917 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5275   -1.4179   -1.5216 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0692   -0.1685    0.1649 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0064   -0.9882    0.8930 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6740    0.7185    0.6768 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5354    0.0937   -0.2054 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6019    0.8495   -0.9960 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2308   -1.1212   -0.6857 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.3496   -1.7465    0.1149 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6518   -1.6134   -1.3602 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2597    0.6597    1.0215 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.3988    1.8397    1.3031 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7152   -0.2977    1.8646 O   0  0  0  0  0  0  0  0  0  0  0  0\n    4.1491    0.2188    2.5749 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  4  1  0\n  2  3  1  0\n  4  5  1  0\n  4  8  2  3\n  5  7  1  0\n  5  6  1  0\n  8  9  1  0\n  9 10  1  0\n  9 12  1  0\n  9 11  1  0\n 12 15  1  0\n 12 13  1  0\n 12 14  1  0\n 15 18  1  0\n 15 17  1  0\n 15 16  1  0\n 18 20  1  0\n 18 23  1  0\n 18 19  1  0\n 20 21  1  0\n 20 22  1  0\n 23 25  1  0\n 23 24  2  0\n 25 26  1  0\nM  END\n","mol");
	viewer_16816845087131624.setStyle({"stick": {}});
	viewer_16816845087131624.zoomTo();
viewer_16816845087131624.render();
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

Visualize(arg_monomer.LHS).infer_bond_orders().create_2D_image('tests/samples/arginine_LHS.png',(400,200))
Visualize(arg_monomer.link).infer_bond_orders().create_2D_image('tests/samples/arginine_link.png',(400,200))
Visualize(arg_monomer.RHS).infer_bond_orders().create_2D_image('tests/samples/arginine_RHS.png',(400,200))
from IPython.display import Image
from IPython.core.display import HTML
arginines = ['arginine_LHS.png','arginine_link.png', 'arginine_RHS.png']
html = '<div style="display:flex">'
for image in arginines:
    html += f'<figure><img src="tests/samples/{image}" style="margin:0 10px" width="300"><figcaption>{image}</figcaption></figure>'
html += '</div>'
display(HTML(html))

```

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

To convert this notebook to markdown, run the following command: 
    ```
    jupyter nbconvert README.ipynb --to markdown
    ```


```python

```
