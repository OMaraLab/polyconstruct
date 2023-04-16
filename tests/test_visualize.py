from polytop.monomer import Monomer
from polytop.topology import Topology
from polytop.visualize import Visualize


def test_2D_monomer():
    arg = Topology.from_ITP('tests/samples/arginine.itp')
    start = arg.get_bond_by_name('N3','H20')
    end = arg.get_bond_by_name('C11','O1')
    arg_monomer = Monomer(arg,start,end)
    arg_monomer.save('tests/samples/arginine_monomer.json')

    Visualize(arg_monomer.LHS).infer_bond_orders().create_2D_image('tests/samples/arginine_LHS.png',(400,200))
    Visualize(arg_monomer.link).infer_bond_orders().create_2D_image('tests/samples/arginine_link.png',(400,200))
    Visualize(arg_monomer.RHS).infer_bond_orders().create_2D_image('tests/samples/arginine_RHS.png',(400,200))
