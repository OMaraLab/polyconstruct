from pathlib import Path
from polytop.junction import Junction, Junctions
from polytop.monomer import Monomer
from polytop.topology import Topology
from polytop.visualize import Visualize
import os
import py3Dmol

def test_visualize_rename_atoms(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir/'arginine.itp')
    Visualize.topology(arg).draw2D(output_dir/'arginine_before_rename.png', show_atom_ID=True, remove_explicit_H=False, size=(400,200))
    arg.auto_rename_atoms()
    Visualize.topology(arg).draw2D(output_dir/'arginine_after_rename.png', show_atom_ID=True, remove_explicit_H=False, size=(400,200))

def test_2D_monomer(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir/'arginine.itp')
    junctions = Junctions()
    arg_monomer = Monomer(arg,[arg.junction('N3','H20').named('N'),arg.junction('N6','H23').named('N'),arg.junction('C11','O1').named('C')])
    arg_monomer.save(output_dir/'arginine_monomer.json')
    
    # delete image first if it exists
    image_path = output_dir/'arginine_monomer.png'
    # if os.path.exists(image_path):
    #     os.remove(image_path)

    Visualize.monomer(arg_monomer).draw2D(image_path,size=(800,300),show_legend=True)

    # assert image was created
    assert os.path.exists(image_path)

