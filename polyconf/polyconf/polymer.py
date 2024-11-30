# TODO: hand polymer ITP connectivity info to polyconf, instead of relying on input pdb's having connectivity
 
#!/usr/bin/env python
from .monomer import Monomer
from .PDB import PDB
 
import numpy as np
import pandas as pd
from tqdm import tqdm

import MDAnalysis as mda
from MDAnalysis.analysis import distances
from math import degrees

import random
import networkx as nx

class Polymer:
    """
    Docstrings go brr
    """
    def __init__(self, firstMonomer: Monomer) -> None:
        """
        Initiate the polymer by building its first monomer. Takes a copy of its
        atoms and sets the residue resids to 1.

        Args:
            firstMonomer (MDAnalysis Universe): first monomer of the polymer
        """
        # Guard clauses to ensure that firstMonomer.residues.resids is a sane call
        if firstMonomer is None:
            raise ValueError('firstMonomer must be a Monomer object')
        if not isinstance(firstMonomer, Monomer):
            raise TypeError('firstMonomer must be a Monomer object')
        
        firstMonomer.residues.resids = 1
        self.first = firstMonomer
        self.polymer = self.first
        self.atoms = self.polymer.atoms

    def select_atoms(self, selection):
        """
        Selection method that selects atoms from the polymer's MDAnalysis
        Universe by leveraging the MDAnalysis atom selection.

        Args:
            selection (str): MDAnalysis atom selection string, for more 
                    information on supported atom selection formats see 
                    [MDAnalysis Atom selection language](https://userguide.mdanalysis.org/stable/selections.html)
        """
        return self.polymer.select_atoms(selection)
    
    def renamer(self, resid, namein, nameout='X'):
        """
        Change selected atom names to a new name like X1, X2, X3. Intended to
        flag dummy atoms for removal. Selected atoms are given a basename, 
        e.g. 'X' defined by the nameout argument, as well as a number. 

        Args:
            resid (int): which residue number to select and rename atoms from
            namein (str): current name of atoms, used to select target dummy atoms
            nameout (str): new basename of atoms selected from resid and with
                    atom name defined by namein
        """
        i = len(self.polymer.select_atoms(f"resid {resid} and name {nameout}*").atoms.names) + 1 # count existing dummy atoms
        dummies=self.polymer.select_atoms(f"resid {resid} and name {namein}")
        for x in dummies.atoms.indices:
            self.polymer.select_atoms(f"resid {resid} and name {namein} and index {x}").atoms.names=[nameout+str(i)]

    def newresid(self):
        """
        Generates a resid that is one larger than the highest existing resid in
        the polymer. Excellent for finding what the number of the next residue
        should be set to, for the nn argument to 'extend()'.

        Returns:
            nn (int): resid one greater than the polymer's current highest resid
        """
        nn = max(self.polymer.residues.resids) + 1
        return nn

    def extend(self, monomer, n, nn, names, joins, ortho=[1,1,1], 
               linearise=False, beta=0): 
        """
        Extend the polymer by adding a monomer

        Args:
            monomer (Monomer): the monomer to add to the polymer
            n (int): residue to extend from
            nn (int): residue to extend to (i.e. the new residue)
            names (dict): four single 'key:value' pairs with keys 
                    P1, Q1, P2 and Q2 for the four dummy atoms that are 
                    'overlaid' during polymer extension. P1 is the real atom in
                    the polymer while P2 is its corresponding dummy atom in the
                    incoming monomer. Similarly, Q1 is the real atom in the 
                    polymer while Q2 is the corresponding dummy atom in the 
                    incoming polymer. Such that P1 and Q1 and bonded, and P2 
                    and Q2 have an equivalent bond. The result will be a new 
                    bond between Q1 and the next atom beside Q2 that is NOT P2.
            joins (list of a 2 item tuple): defines which two atoms of residue
                    'n' (first value) and residue 'nn' (second value) will 
                    be bonded
            ortho (list with 3 values): used for linearise and defines which 
                    axis to align along (e.g. align along x with ortho=[1,0,0])
                    Default value is '[1,1,1]'.
            linearise (bool): make the polymer linear. Good for debugging and
                    visually confirming structure is correct. 
                    Default value is False.
            beta (int or float): beta is used to group monomers into 
                    categories, so that each branch can be separated. Beta is
                    analogous to segments. The tempfactors attribute of all of
                    the atoms in the incoming monomer are set to the value of 
                    beta. Default value is 0.
            
        
        For example: {'P1':'CA','Q1':'C','P2':'CMA','Q2':'CN'}, with n=1 and 
                nn=2 and joins=[('C','CA')] would produce a polymer of 
                '...-CA1-C1-CA2-C2-...'

        Further notes:
            * n and nn enable branching by specifying from and to connections between monomers and beta is analogous to segments
            * extend a polymer u by a monomer u_, by fitting the backbone atoms (P2, Q2) from the new monomer to (P1,Q1) from the existing residue n, 
                then joins the monomer nn to the existing universe with a bond between each pair of atoms in joins
        
            ATOM NOMENCLATURE:
            * P1 and Q1 are atoms in monomer n
            * P2 and Q2 are dummy atoms in monomer nn, which correspond to the atoms P1 and Q1
            * R and S are atoms in monomer nn that describe a vector to be aligned to the ortho vector during the linearise fit; 
                intended for making straight orthoganal branches for ease of visibility 
                this must be chosen carefully so the rotation doe not break stereochemistry or alignment
            JOINS NOMENCLATURE:
            * joins contains pairs of atoms to be linked, of the form (X_n,Y_nn)
            * X_n   is some atom in residue n, the final residue before extension
            * Y_n+1 is some atom in residue nn, the new residue
        
        Note: this function preserves all dummy atoms. Removing dummy atoms 
            during extension causes substantial problems with indexing. They 
            must only be removed once the polymer is fully built. Save your 
            polymer to a .gro file with the PDB.save() method, which will strip
            any dummy atoms.
        """

        P1 = self.polymer.select_atoms(f'resid {n} and name {names["P1"]}').positions[-1]
        Q1 = self.polymer.select_atoms(f'resid {n} and name {names["Q1"]}').positions[-1]

        u_ = monomer # monomer = mdict['path'][monomer]
        u_.atoms.tempfactors=float(beta) 

        u_.residues.resids = nn

        P2 = u_.select_atoms(f'resid {nn} and name {names["P2"]}').positions[0]

        # first, move CMA_n+1 to CA_N
        T = P1 - P2

        u_.atoms.translate(T)

        # next, rotate around cross product of backbone vectors to align C_n to CN_n+1
        P2 = u_.select_atoms(f'resid {nn} and name {names["P2"]}').positions[0]
        Q2 = u_.select_atoms(f'resid {nn} and name {names["Q2"]}').positions[0]

        v1 = Q2 - P2
        v1_n = np.linalg.norm(v1)

        v2 =  Q1 - P1
        v2_n = np.linalg.norm(v2)

        theta = degrees(np.arccos(np.dot(v1,v2)/(v1_n * v2_n)))

        k = np.cross(v1,v2)
        u_r1 = u_.atoms.rotateby(theta,axis=k,point=P1)

        if linearise:
            R= u_.select_atoms('resid '+str(nn)+' and name '+names['R']).positions[0]
            S= u_.select_atoms('resid '+str(nn)+' and name '+names['S']).positions[0]
            RS=S-R
            RS_n=np.linalg.norm(RS)
            ortho_n=np.linalg.norm(ortho)
            theta = degrees(np.arccos(np.dot(RS,ortho)/(RS_n * ortho_n)))

            k=np.cross(RS,ortho)
            
            u_r2 = u_r1.atoms.rotateby(theta,axis=k,point=R)

            R= u_.select_atoms('resid '+str(nn)+' and name '+names['R']).positions[0]
            S= u_.select_atoms('resid '+str(nn)+' and name '+names['S']).positions[0]
            RS=S-R
            RS_n=np.linalg.norm(RS)
            ortho_n=np.linalg.norm(ortho)
            check = degrees(np.arccos(np.dot(RS,ortho)/(RS_n * ortho_n)))
            if abs(check)>1: print('linear check =' ,check)

            u_r1=u_r2

        # combine extended polymer into new universe
        new = mda.Merge(self.polymer.atoms, u_r1.atoms)

        # add new bonds linking pairs of atoms (X_n,Y_nn) 
        for pair in joins:
            X = new.select_atoms("resid "+str(n)+" and name "+pair[0]).indices[0]
            Y = new.select_atoms("resid "+str(nn)+" and name "+pair[1]).indices[0]
            new.add_bonds([(X,Y)])

        new.dimensions = list(new.atoms.positions.max(axis=0) + [0.5,0.5,0.5]) + [90]*3
        self.polymer = new.copy()
    
    #TODO: REPLACE WITH NAC GENCONF VERSION??
    def genconf(self, atomPairs, dummies, length, runs = 5, attempts = 20, cutoff = 0.8, fname = 'polymer_conf', pdb = False):
        """
        Rotate bonds in the polymer to generate various conformations without atom overlap.

        Args:
            atomPairs (list of 2 item tuples): pairs of bonded atoms to rotate
                    around, only 1-2 pairs of bonded atom types in the polymer 
                    should be provided.
            dummies (list): list of names of dummy atoms
            length (int): length of the polymer (i.e. the number of monomers)
            runs (int): number of conformations to generate, default is 5
            attempts (int): maximum number of shuffle rotations to attempt with
                    clashes being detected before retrying the shuffle from the
                    start, default is 20
            cutoff (float): minimum distance allowed between any two atoms, the
                    default is 0.8
            fname (str): base file name, default is 'polymer_conf'
            pbd (bool): save output as PDB if True, else save as default GROMACS

        Writes each of the generated comformations to their own numbered output file
        """
        j=1
        conf = self.polymer.copy()
        while j <= runs:
            tries = 0
            for i, index in enumerate(self.polymer.residues, start=1):     
                for a1,a2 in atomPairs:
                    if a1[-1]=='+' or a2[-1]=='+':
                        sign = '+'
                    elif a1[-1]=='-' or a2[-1]=='-':
                        sign = '-'
                    else:
                        sign = ''

                    if (not sign): # [CA C]
                        conf, a, b = self.shuffle(conf,a1, i,a2, i)
                    elif (sign == '+'): # [CA C+]
                        if (i != length):
                            conf, a, b =self.shuffle(conf,a1, i,a2, i+1)
                    else: # [CA C-]
                        if (i != 1):
                            conf, a, b =self.shuffle(conf,a1, i,a2, i-1)

                    #convert a and b from atom indicies to atomgroups
                    aAtoms = conf.atoms[a.nodes()]
                    bAtoms = conf.atoms[b.nodes()]
                    dist=distances.distance_array(aAtoms.atoms.positions, bAtoms.atoms.positions).min() # this is the clash detection
                    if dist >= cutoff:
                        print(f'genconf {i} success')
                        toSave = PDB(conf)
                        toSave.cleanup()
                        toSave.save(dummies, fname = f'{fname}_{i}', selectionString=None, pdb=pdb)
                        j+=1
                        tries = 0
                        #proceed to next - save at end of a run only
                    else:
                        #redo shuffle (until max trieds reached, then give up)
                        if (tries == attempts):
                            print('max number of attempts reached; restarting shuffle')
                            tries = 0
                            continue
                        else:
                            print(f'genconf {i} clash; retrying')
                            tries+=1
   
    def shuffle(self, u, a1, a1resid, a2, a2resid, mult=3):
        """
        based on a tutorial by richard j gowers; http://www.richardjgowers.com/2017/08/14/rotating.html
        """ 
        pair = self.polymer.select_atoms(f'(resid {a1resid} and name {a1}) or (resid {a2resid} and name {a2})' ) 
        bond = self.polymer.atoms.bonds.atomgroup_intersection(pair,strict=True)[0]
        g = nx.Graph()
        g.add_edges_from(self.polymer.atoms.bonds.to_indices()) # get entire residue as a graph
        g.remove_edge(*bond.indices)     # remove the bond from the graph
        # unpack the two unconnected graphs
        a, b = (nx.subgraph(g,c) for c in nx.connected_components(g))
        # call the graph without the CA atom the 'head', this will be rotated
        fore_nodes = a
        fore = self.polymer.atoms[fore_nodes.nodes()]
        v = bond[1].position - bond[0].position
        o = (fore & bond.atoms)[0]
        rot = random.randrange(0,mult)*int(360/mult) # rotate fore by a random multiplicity
        fore.rotateby(rot, v, point=o.position)
        return(u, a, b)
    
    
    def gencomp(mdict, length, fill, middle, count = True, frac = False): 
        """
        Generate the linear conformation of a polymer from a given length, 
        specified monomer composition and a dictionary of monomers.

        Args:
            mdict (dict): dictionary of resnames in the middle with how many
                    of each type
            length (int): desired length of the polymer 
                    (i.e. number of monomers)
            fill (list): list of monomers which can be used to build the 
                    polymer middle (i.e. are not terminal monomers) if there 
                    are not enough in 'middle' to complete the monomer
            middle (list): list of monomers to build the polymer from
            count (bool): generate absolute composition with counts (e.g. 4 of
                    monomer A and 12 of monomer B), default is True
            frac (bool): generate relative composition with fractions (e.g. 25%
                    monomer A and 75% monomer B), default is False
        
        Returns:
            polyList: a list with the sequence of randomly ordered monomers of
                    the desired length and composition to build the polymer
        """
        #TODO: is this able to generate instructions for branched polymers?
        polycomp = []

        if count:
            for m in middle:
                rcount = mdict['count'][m]
                for i in range(rcount):
                    polycomp += [m] # build a list with one of each monomer unit present in the final polymer

        elif frac:
            for m in middle:
                rfrac = int(length * mdict['frac'][m]) # fraction of total length that is monomer x
                for i in range(rfrac):
                    polycomp += [m] # build a list with one of each monomer unit present in the final polymer
        
        # then randomise positions of monomers
        # 'polycomp' is ordered list of X number of monomers -> will add all of monomer A required, before doing all of B, etc.

        filler = [x for x in fill]
        while len(polycomp) < length:
            sel = random.sample(filler,1)
            polycomp += sel # if the polymer is shorter than the desired length, continue to fill from the given list 
            filler = [x for x in filler if not x in sel]
            if len(filler) == 0: 
                filler = [x for x in fill]

        polyList=random.sample(polycomp,length) # reorder the list of monomers randomly

        polyList[0] = polyList[0][:-1] + 'I' # convert first monomer from middle to initiator -> add extra bit onto middle mono
        # TODO: REWORK THIS TO BUILD IN SAME WAY AS POLYTOP!

        polyList[-1] = polyList[-1][:-1] + 'T' # convert last monomer from middle to terminator -> add extra bit onto middle mono
        return(polyList) # list of length items with all monomers in order -> used to build w PDBs
