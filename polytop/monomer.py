from functools import singledispatchmethod

class Monomer():
    def __init__(self,molecule,LHS,RHS) -> None:
        self.molecule= molecule
        self.LHS = LHS
        self.RHS = RHS
        pass

