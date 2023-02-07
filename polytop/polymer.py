from functools import singledispatchmethod
from enum import Enum

class distributePartialCharges(Enum):
    ACROSSMONOMER = 1
    ACROSSPOLYMER = 2
    NEARESTMOSTCHARGEDATOM = 3
    
class Polymer():
    def __init__(self,distributeCharges: distributePartialCharges ) -> None:
        self.distCharges = distributeCharges
        return None

