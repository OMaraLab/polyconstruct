from pathlib import Path
from polytop.junction import Junction, Junctions
from polytop.monomer import Monomer
from polytop.topology import Topology
from polytop.visualize import Visualize
import os



def test_2D_monomer(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir/'arginine.itp')
    junctions = Junctions()
    junctions.add(Junction('N',arg.get_bond('N3','H20')))
    junctions.add(Junction('N',arg.get_bond('N6','H23')))
    junctions.add(Junction('C',arg.get_bond('C11','O1')))
    arg_monomer = Monomer(arg,junctions)
    arg_monomer.save(output_dir/'arginine_monomer.json')
    
    # delete image first if it exists
    image_path = output_dir/'arginine_monomer.png'
    if os.path.exists(image_path):
        os.remove(image_path)

    Visualize.monomer(arg_monomer).draw2D(image_path,(400,200))

    # assert image was created
    assert os.path.exists(image_path)
