from functools import singledispatchmethod

class Molecule():
    def __init__(self) -> None:
        self.atoms =[]
        self.pairs = []
        self.triplets=[]
        return None

    @classmethod
    def fromITP(cls, filename: str):
        newMolecule = cls()
        # read text file
        return newMolecule
