import datetime
import os
import warnings


class comment:
    def __init__(self):
        self.content = None
        self.isParam = False
        pass

    def parse(self, line):
        self.content = line.rstrip("\n")
        if len(self.content) > 1:
            self.isParam = self.content[1] == "@"
        return self

    def summary(self):
        return " {}\n".format(self.content)


class molecule:
    def __init__(self):
        self.molname = None
        self.nrexcl = None
        self.mapping = None
        pass

    def parse(self, line):
        columns = line.split()
        self.molname = columns[0]
        self.mapping = columns[0]
        self.nrexcl = columns[1]
        return self

    def summary(self):
        return " {}\t{}\n".format(self.molname, self.nrexcl)


class atom:
    def __init__(self):
        self.id = None
        self.type = None
        self.resnr = None
        self.residu = None
        self.atom = None
        self.cgnr = None
        self.charge = None
        self.comment = None
        pass

    def parse(self, line):
        columns = line.split()
        self.id = int(columns[0])
        self.type = columns[1]
        self.resnr = columns[2]
        self.residu = columns[3]
        self.atom = columns[4]
        self.cgnr = columns[5]
        self.charge = columns[6]
        if ";" in line:
            self.comment = line.split(";")[-1]
        if not any(chr.isdigit() for chr in columns[4]):
            warnings.warn(f"No index in atom {columns[4]}, atom id {self.id}. Please check your ITP file.")
        return self

    def summary(self):
        if not self.comment:
            return " {}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                self.id,
                self.type,
                self.resnr,
                self.residu,
                self.atom,
                self.cgnr,
                self.charge,
            )
        else:
            return " {}\t{}\t{}\t{}\t{}\t{}\t{} ; {}\n".format(
                self.id,
                self.type,
                self.resnr,
                self.residu,
                self.atom,
                self.cgnr,
                self.charge,
                self.comment,
            )

    def __repr__(self):
        return self.atom

    def __str__(self):
        return self.atom


class bond:
    def __init__(self):
        self.i = None
        self.j = None
        self.funct = None
        self.length = None
        self.force = None
        self.comment = None
        pass

    def parse(self, line):
        columns = line.split()
        self.i = int(columns[0])
        self.j = int(columns[1])
        self.funct = columns[2]
        self.length = columns[3]
        self.force = columns[4]
        if ";" in line:
            self.comment = line.split(";")[-1]
        return self

    def summary(self):
        if not self.comment:
            return "{}\t{}\t{}\t{}\t{}\n".format(
                self.i, self.j, self.funct, self.length, self.force
            )
        else:
            return "{}\t{}\t{}\t{}\t{}\t ; {}\n".format(
                self.i, self.j, self.funct, self.length, self.force, self.comment
            )

    @classmethod
    def create(cls, i, j, funct, length, force, comment=None):
        newbond = cls()
        newbond.i = i
        newbond.j = j
        newbond.funct = funct
        newbond.length = length
        newbond.force = force
        newbond.comment = comment
        return newbond


class angle:
    def __init__(self):
        self.i = None
        self.j = None
        self.k = None
        self.funct = None
        self.angle = None
        self.force = None
        self.comment = None
        pass

    def parse(self, line):
        columns = line.split()
        self.i = int(columns[0])
        self.j = int(columns[1])
        self.k = int(columns[2])
        self.funct = columns[3]
        self.angle = columns[4]
        self.force = columns[5]
        if ";" in line:
            self.comment = line.split(";")[-1]
        return self

    def summary(self):
        if not self.comment:
            return "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                self.i, self.j, self.k, self.funct, self.angle, self.force
            )
        else:
            return "{}\t{}\t{}\t{}\t{}\t{}\t ; {}\n".format(
                self.i, self.j, self.k, self.funct, self.angle, self.force, self.comment
            )

    @classmethod
    def create(cls, i, j, k, funct, angle, force, comment=None):
        newangle = cls()
        newangle.i = i
        newangle.j = j
        newangle.k = k
        newangle.funct = funct
        newangle.angle = angle
        newangle.force = force
        newangle.comment = comment
        return newangle


class ITP:
    def __init__(self):
        self.molecule = None
        self.comments = []
        self.atoms = []
        self.bonds = []
        self.angles = []
        pass

    def __str__(self):
        if self.molecule == None:
            return "empty ITP"
        else:
            return "ITP for molecule:{} - {} atoms".format(
                self.molecule.molname, len(self.atoms)
            )

    def load(self, path):
        tag = None
        with open(path, "r") as input:
            for line in input.readlines():
                if not line.strip():
                    continue
                if line[0] == ";":
                    if tag == None:
                        self.comments.append(comment().parse(line))
                    continue
                if line[0] == "[":
                    # remove spaces from line
                    tag = line.replace(" ", "")
                    continue
                if tag == "[moleculetype]\n":
                    self.molecule = molecule().parse(line)
                elif tag == "[atoms]\n":
                    self.atoms.append(atom().parse(line))
                elif tag == "[bonds]\n":
                    self.bonds.append(bond().parse(line))
                elif tag == "[angles]\n":
                    self.angles.append(angle().parse(line))
        return self

    def status(self):
        print("molecule:{}".format(self.molecule.molname))
        print("{} atoms".format(len(self.atoms)))
        print("{} bonds".format(len(self.bonds)))
        print("{} angles".format(len(self.angles)))
        return self

    def remap(self, newID):
        for atom in self.atoms:
            atom.id += newID
        for bond in self.bonds:
            bond.i += newID
            bond.j += newID
        for angle in self.angles:
            angle.i += newID
            angle.j += newID
            angle.k += newID
        return self

    def appendITP(self, path, name, mapping, charge):
        # find HWM
        hwm = 0
        for atom in self.atoms:
            if atom.id > hwm:
                hwm = atom.id
        newITP = ITP().load(path)
        newITP.remap(hwm)
        self.atoms.extend(newITP.atoms)
        self.atoms[0].type = "Nda"
        self.atoms[0].comment = "linking bead, note Nda"
        self.atoms[1].comment = "first phosphatidyl group"
        self.atoms[2].comment = "first glycerol group"
        self.atoms[hwm].comment = "second phosphatidyl group"
        self.atoms[hwm + 1].comment = "second glycerol group"
        for atom in self.atoms:
            if atom.id > hwm:
                atom.atom = atom.atom.replace("A", "C")
                atom.atom = atom.atom.replace("B", "D")
                if atom.atom == "PO4":
                    atom.atom = "PO42"
                elif atom.atom == "GL1":
                    atom.atom = "GL3"
                elif atom.atom == "GL2":
                    atom.atom = "GL4"
        self.atoms[0].atom = "GL5"
        self.atoms[1].atom = "PO41"

        self.bonds.extend(newITP.bonds)
        self.angles.extend(newITP.angles)

        self.bonds[0].length = 0.37
        self.bonds[0].force = 5500
        self.bonds.append(
            bond.create(
                1, hwm + 1, funct=1, length=0.37, force=5500, comment="UA 0.372 (avg)"
            )
        )
        self.angles.append(
            angle.create(
                2,
                1,
                hwm + 1,
                funct=2,
                angle=105,
                force=45,
                comment="linking bead: PO41-GL5-PO42",
            )
        )
        self.angles.append(
            angle.create(
                1,
                hwm + 1,
                hwm + 2,
                funct=2,
                angle=110,
                force=25,
                comment="linking bead to second glycerol group: PO41-GL5-PO42",
            )
        )
        self.angles.append(
            angle.create(
                1,
                2,
                3,
                funct=2,
                angle=110,
                force=25,
                comment="linking bead to first glycerol group:  PO41-GL5-PO42",
            )
        )
        newcomments = []
        newcomments.append(
            ";;;; Martini lipid topology for Cardiolipin {}, mapped to code {}, generated using \n".format(
                name, mapping
            )
        )
        newcomments.append("\n")
        newcomments.append("\n")
        newcomments.append("; Description:\n")
        newcomments.append(
            ";   A general model Cardiolipin using the following beads \n"
        )
        atomdesc = ""
        for atom in self.atoms:
            atomdesc += atom.atom + " "
        newcomments.append("; " + atomdesc + "\n")
        newcomments.append("\n")
        newcomments.append("; Parameterization:\n")
        newcomments.append(
            ";   This topology follows the standard Martini 2.0 lipid definitions and building block rules.\n"
        )
        newcomments.append("\n")
        newcomments.append("; Reference(s): \n")
        newcomments.append(
            ";   S.J. Marrink, A.H. de Vries, A.E. Mark. Coarse grained model for semi-quantitative lipid simulations. JPC-B, 108:750-760, \n"
        )
        newcomments.append(";   2004. doi:10.1021/jp036508g \n")
        newcomments.append(
            ";   S.J. Marrink, H.J. Risselada, S. Yefimov, D.P. Tieleman, A.H. de Vries. The MARTINI force field: coarse grained model for \n"
        )
        newcomments.append(
            ";   biomolecular simulations. JPC-B, 111:7812-7824, 2007. doi:10.1021/jp071097f \n"
        )
        newcomments.append(
            ";   T.A. Wassenaar, H.I. Ingolfsson, R.A. Bockmann, D.P. Tieleman, S.J. Marrink. Computational lipidomics with insane: a versatile \n"
        )
        newcomments.append(
            ";   tool for generating custom membranes for molecular simulations. JCTC, 150410125128004, 2015. doi:10.1021/acs.jctc.5b00209\n"
        )
        newcomments.append(
            ";   M. Dahlberg. Polymorphic phase behavior of cardiolipin derivatives studied by coarse-grained molecular dynamics. JPC-B, \n"
        )
        newcomments.append(";   111:7194â€“7200, 2007. doi:10.1021/jp071954f\n")
        newcomments.append(
            ";   M. Dahlberg, A. Maliniak. Mechanical properties of coarse-grained bilayers formed by cardiolipin and zwitterionic lipids.\n"
        )
        newcomments.append(";   JCTC, 6:1638-1649, 2010. doi:10.1021/ct900654e\n")

        newcomments.append("\n")
        newcomments.append(
            "; Created by linking {} generated with the following parameters \n".format(
                self.molecule.molname
            )
        )
        for item in self.comments:
            if item.isParam:
                newcomments.append(item.content)
        newcomments.append("\n")
        newcomments.append(
            "; with {} generated with the following parameters \n".format(
                newITP.molecule.molname
            )
        )
        for item in newITP.comments:
            if item.isParam:
                newcomments.append(item.content)
        newcomments.append("\n")
        newcomments.append(
            "; Generated using the Martini lipid itp generator version 0.6\n"
        )
        newcomments.append("\n")
        newcomments.append("; Generated on {}\n".format(datetime.datetime.now()))
        self.comments.clear()
        for newcomment in newcomments:
            self.comments.append(comment().parse(newcomment))
        self.molecule.molname = name
        self.molecule.mapping = mapping
        for atom in self.atoms:
            atom.residu = mapping
            atom.cgnr = atom.id
        if charge == -2:
            self.atoms[1].charge = -1
            self.atoms[hwm].charge = -1
        elif charge == -1:
            self.atoms[1].charge = -1
            self.atoms[hwm].charge = 0
        else:
            self.atoms[1].charge = 0
            self.atoms[hwm].charge = 0
        return self

    def write(self, path):
        with open(path, "w") as output:
            for comment in self.comments:
                output.write(comment.summary())
            output.write("[moleculetype]\n")
            output.write("; molname      nrexcl\n")
            output.write("{}\t{}\n".format(self.molecule.mapping, self.molecule.nrexcl))
            output.write("\n[atoms]\n")
            output.write("; id	type 	resnr 	residu 	atom 	cgnr 	charge\n")
            for atom in self.atoms:
                output.write("  {}".format(atom.summary()))
            output.write("\n[bonds]\n")
            output.write("; i	j 	funct 	length 	force.c.\n")
            for bond in self.bonds:
                output.write("  {}".format(bond.summary()))
            output.write("\n[angles]\n")
            output.write("; i	j	k 	funct 	angle 	force.c.\n")
            for angle in self.angles:
                output.write("  {}".format(angle.summary()))
        return self
