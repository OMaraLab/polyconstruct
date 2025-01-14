from typing import List

class Junction:
    """
    Junctions are the polymerization sites of a monomer topology. They are
    defined by a name and a list of Bonds.
    """
    from .Atoms import Atom
    from .Bonds import Bond

    #TODO: make name enforcably not optional
    def __init__(self, monomer_atom: Atom, residue_atom: Atom, name: str = None):
        """
        Junctions are the polymerization sites of a monomer topology. They are
        defined by a name and a list of Bonds. The name should be unique to the
        Junction such that it can be unambigously identified by its name, thus
        ensuring that polymers can be built with consistent, desired
        conformations.

        :param monomer_atom: the Atom which will remain with this Monomer after
                polymerisation, and will obtain a new bond.
        :type monomer_atom: Atom
        :param residue_atom: the Atom which will be lost during polymerisation,
                analogous to the leaving atom in a chemical reaction.
        :type residue_atom: Atom
        :param name: the unique name of the Junction type, defaults to None.
        :type name: str, optional
        :raises ValueError: if a Junction cannot be determined from the two
                provided atoms (i.e. if they are not bonded).
        """
        if name is None:
            name = f"{monomer_atom.atom_name}-{residue_atom.atom_name}"
        self.name = name
        self.monomer_atom = monomer_atom
        self.residue_atom = residue_atom
        if not residue_atom in monomer_atom.bond_neighbours():
            raise ValueError("Junction location cannot be found")
        
    def named(self, newname: str) -> "Junction":
        """
        Use to name or rename a Junction after it has been created.

        :param newname: the desired, unique name that will be used to identify
                this Junction.
        :type newname: str
        :return: the Junction object with a name
        :rtype: Junction

        Depreciation warning:
        * This function will be depreciated shortly in favour of enforced
            setting the name attribute of a Junction when it is created.

        Instead of:   junction = Junction(atomA, atomB).named("name")
        Use:          junction = Junction(atomA, atomB, name="name")
        """
        self.name = newname
        return self
        
    def to_dict(self) -> dict:
        """
        Convert this Junction to a dictionary representation.

        The structure of the dictionary is as below:
        {"name": self.name,
        "monomer_atom": self.monomer_atom.atom_name,
        "residue_atom": self.residue_atom.atom_name}

        :return: a dictionary containing the names of its Atoms and the unique
                name of this Junction.
        :rtype: dict
        """
        return {
            "name": self.name,
            "monomer_atom": self.monomer_atom.atom_name,
            "residue_atom": self.residue_atom.atom_name
        }

    @classmethod
    def from_dict(cls, data: dict, atoms: List[Atom]):
        """
        Create a new Junction from a dictionary, such as that created with
        Junction.to_dict().

        The structure of the dictionary is as below:
        {"name": self.name,
        "monomer_atom": self.monomer_atom.atom_name,
        "residue_atom": self.residue_atom.atom_name}

        :param data: dictionary containing data to make a Junction, generate
                with 'to_dict()'.
        :type data: dict
        :param atoms: list of Atoms. The list may contain more than 2 atoms, as
                long as the two atoms specified in the data dict are present.
        :type atoms: List[Atom]
        :return: a new Junction
        :rtype: Junction
        """
        name = data["name"]
        from .Atoms import Atom
        monomer_atom_name = data["monomer_atom"]
        monomer_atom = next(atom for atom in atoms if atom.atom_name == monomer_atom_name)
        residue_atom_name = data["residue_atom"]
        residue_atom = next(atom for atom in atoms if atom.atom_name == residue_atom_name)
        return cls(monomer_atom,residue_atom,name)

    @classmethod
    def from_topology(cls, topology: "Topology", monomer_atom_name: str, 
                      residue_atom_name: str, residue_id: int = None, 
                      name: str = None) -> "Junction":
        """
        Create a new Junction from a monomer or polymer Topology, and the names
        and shared residue id of the two Atoms used to form the Junction.

        :param topology: a polymer or monomer Topology. Obtain from the
                '.topology' attribute of either a Monomer or Polymer.
        :type topology: Topology
        :param monomer_atom_name: the name of the Atom which will remain with
                this Topology after polymerisation, and will obtain a new bond.
        :type monomer_atom_name: str
        :param residue_atom_name: the name of the Atom which will be lost
                during polymerisation, analogous to the leaving atom in a
                chemical reaction.
        :type residue_atom_name: str
        :param residue_id: the id number of the residue both atoms belong to,
                defaults to None
        :type residue_id: int, optional
        :param name: the unique name used to identify the new Junction,
                defaults to None
        :type name: str, optional
        :return: a new Junction
        :rtype: Junction
        """
        monomer_atom = topology.get_atom(monomer_atom_name, residue_id)
        residue_atom = topology.get_atom(residue_atom_name, residue_id)
        return cls(monomer_atom, residue_atom, name)

    def __repr__(self) -> str:
        return f"(\"{self.name}\":{self.monomer_atom.atom_name}-{self.residue_atom.atom_name})"
    
class Junctions(list):
    """
    Class used to group and list Junctions of a Monomer or Polymer for RDKit
    visualisation. Used by the Visualize Class.
    """
    def add(self, junction: Junction):
        """
        Add a new Junction to the Junction list.

        :param junction: a Junction to add to the Junction list.
        :type junction: Junction
        """
        self.append(junction)

    def get_junctions(self) -> "Junctions":
        """
        Getter method to access the list of Junctions.

        :return: the list of Junction objects (i.e. a 'Junctions' object).
        :rtype: Junctions
        """
        return self

    def remove(self, junction: Junction):
        """
        Remove a Junction from the list of Junctions.

        :param junction: the Junction object to be removed from the Junctions.
        :type junction: Junction
        """
        super().remove(junction)

    def named(self, name: str) -> list[Junction]:
        """
        Returns a list of Junctions with the name provided.

        :param name: the name of a Junction/s to retrieve.
        :type name: str
        :return: list[Junction]
        :rtype: _type_
        """
        return [junction for junction in self if junction.name == name]

    def to_dict(self) -> list[dict]:
        """
        Return a list containing the dictionary representations of all of the
        Junctions in the Junctions list.

        The structure of each dictionary within the list is as below:
        {"name": self.name,
        "monomer_atom": self.monomer_atom.atom_name,
        "residue_atom": self.residue_atom.atom_name}

        :return: a list of Junction dictionary representations.
        :rtype: list[dict]
        """
        return [junction.to_dict() for junction in self]

    from .Atoms import Atom
    @classmethod
    def from_dict(cls, data: list, atoms: List[Atom]) -> "Junctions":
        """
        Create a new list of Junctions from a list of Junction dictionary
        representations, such as that created with Junctions.to_dict().

        The structure of each dictionary within the list is as below:
        {"name": self.name,
        "monomer_atom": self.monomer_atom.atom_name,
        "residue_atom": self.residue_atom.atom_name}

        :param data: list of dictionaries containing data to make Junctions,
                generate with 'to_dict()'.
        :type data: list
        :param atoms: list of Atoms. The list may contain more than 2 atoms, as
                long as the two atoms specified in the data dict are present.
        :type atoms: List[Atom]
        :return: a new Junctions (i.e. list of Junction objects)
        :rtype: Junctions
        """
        junctions = cls()
        for junction_dict in data:
            junction = Junction.from_dict(junction_dict, atoms)
            junctions.add(junction)
        return junctions

    def __repr__(self) -> str:
        return f"({len(self)}) {','.join(j.__repr__() for j in self)}"