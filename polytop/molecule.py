from functools import singledispatchmethod

from .ITP import ITP


class Molecule:
    def __init__(self) -> None:
        self.atoms = []
        self.pairs = []
        self.triplets = []
        return None

    @singledispatchmethod
    @classmethod
    def fromITP(cls, arg):
        raise NotImplemented("fromITP only accepts ITP or str(filename) as arguments")

    @fromITP.register
    @classmethod
    def _(cls, ITPObject: ITP):
        newMolecule = cls()
        newMolecule.atoms = ITPObject.atoms
        newMolecule.pairs = ITPObject.bonds
        newMolecule.triplets = ITPObject.angles
        return newMolecule

    @fromITP.register
    @classmethod
    def _(cls, filename: str):
        ITP_file = ITP()
        ITP_file.load(filename)
        return cls.fromITP(ITP_file)

    def set_start_atoms():
        pass

    def set_stop_atoms():
        pass
