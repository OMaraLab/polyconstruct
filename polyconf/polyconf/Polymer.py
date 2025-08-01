#!/usr/bin/env python
from __future__ import annotations
from .Monomer import Monomer
 
import numpy as np
from tqdm import tqdm

import MDAnalysis as mda
from MDAnalysis.analysis import distances
from math import degrees

import random
import networkx as nx

# TODO: hand polymer ITP connectivity info to polyconf, instead of relying on input pdb's having connectivity
# TODO set a random seed

class Polymer:
    """
    Create a Polymer from a Monomer and 'extend' it with more Monomers to
    create a Polymer with the desired structure and connectivity.
    """
    def __init__(self, firstMonomer: Monomer, keep_resids=False) -> None:
        """
        Initiate the polymer by building its first monomer. Takes a copy of its
        atoms and sets the residue resids to 1.

        :param firstMonomer: first monomer of the polymer
        :type firstMonomer: Monomer
        :param keep_resids: if True, do not set the resid of all atoms in
                firstMonomer to 1, defaults to False
        :type keep_resids: bool, optional
        :raises ValueError: if firstMonomer is None
        :raises TypeError: if firstMonomer is not a Monomer object
        """

        # Guard clauses to ensure that firstMonomer.residues.resids is a sane call
        if firstMonomer is None:
            raise ValueError('firstMonomer must be a Monomer object')
        if not isinstance(firstMonomer, Monomer):
            raise TypeError('firstMonomer must be a Monomer object')
        
        if not keep_resids:
            firstMonomer.residues.resids = 1
        self.first = firstMonomer
        self.polymer = self.first
        self.atoms = self.polymer.atoms
        self.dimensions = self.polymer.dimensions

    def select_atoms(self, selection) -> mda.AtomGroup:
        """
        Selection method that selects atoms from the polymer's MDAnalysis
        Universe by leveraging the MDAnalysis atom selection.

        :param selection: MDAnalysis atom selection string, for more 
                    information on supported atom selection formats see 
                    `MDAnalysis Atom selection language <https://userguide.mdanalysis.org/stable/selections.html>`
        :type selection: str
        :return: An MDAnalysis AtomGroup containing only the selected atoms
        :rtype: mda.AtomGroup
        """
        return self.polymer.select_atoms(selection)
    
    def copy(self) -> Polymer:
        """
        Use MDAnalysis Universe.copy() to create a new MDAnalysis Universe that
        is an exact copy of this one containing the polymer. Deepcopy is not
        used as it can be problematic with open file sockets.

        :return: A new Polymer that is an exact copy of this polymer
        :rtype: Polymer
        """
        new = self.polymer.copy()
        return Polymer(Monomer.monomer_from_u(new), keep_resids=True)
    
    def renamer(self, resid, namein, nameout='X'):
        """
        Change selected atom names to a new name.
        Intended to flag dummy atoms for removal. 
        Selected atoms are given a basename defined by the nameout argument, 
        e.g. 'X' , and a number. 

        :param resid: which residue number to select and rename atoms from
        :type resid: int
        :param namein: current name of atoms, used to select target dummy atoms
        :type namein: str
        :param nameout: new basename of atoms selected from resid and with
                    atom name defined by namein, defaults to 'X'
        :type nameout: str, optional
        """
        i = len(self.polymer.select_atoms(f"resid {resid} and name {nameout}*").atoms.names) + 1 # count existing dummy atoms
        dummies=self.polymer.select_atoms(f"resid {resid} and name {namein}")
        for x in dummies.atoms.indices:
            self.polymer.select_atoms(f"resid {resid} and name {namein} and index {x}").atoms.names=[nameout+str(i)]


    def newresid(self):
        """
        Returns maximum resid in the polymer plus 1, thus the next available
        resid number an incoming monomer can be assigned.
        
        :return: the polymer's current highest resid plus one.
        :rtype: int
        """
        n = int(max(self.polymer.residues.resids) + 1)
        return (n)


    def maxresid(self):
        """
        Returns maximum resid in the polymer
        
        :return: the polymer's current highest resid
        :rtype: int
        """
        n = int(max(self.polymer.residues.resids))
        return (n)


    def extend(self, monomer, n, nn, names, joins, ortho=[1,1,1], 
               linearise=False, beta=0):
        """
        Extend the polymer by adding a monomer

        :param monomer: the monomer to add to the polymer
        :type monomer: Monomer
        :param n: residue (identified by its resid) to extend from
        :type n: int
        :param nn: residue (identified by its resid) to extend to 
                (i.e. the new residue)
        :type nn: int
        :param names: a dictionary of single 'key:value' pairs with keys
                    P, Q, R and S for the atoms that are used to define the mapping 
                    during polymer extension.  
                    
                    There are two additional optional values, V1 and V2, which are only required in linear extend.
                    
                    The definition is as follows:
                    
                    * Atoms P and R are a pair of bonded atoms located in the incoming monomer.  
                    * Atoms Q and S are a pair of bonded atoms located in residue n of the existing polymer.
                    * Atoms P and Q are equivalent in the final polymer, and either P or Q is a dummy atom.
                    * Atoms R and S are equivalent in the final polymer, and either R or S is a dummy atom.
                    * Bonds PR and QS are equivalent in the final polymer.
                    `extend()` is agnostic as to which atom in each pair is real, and which atom is a dummy
                    after extension, you can use renamer() to explicilty set dummy atoms
                    * Atoms V1 and V2 a pair of atoms located in the incoming monomer.  
        :type names: dict
        :param joins: defines a pair of atoms to connected after extension
                    the first atom in a pair belongs to residue 'n', and the second belongs to residue 'nn'
        :type joins: list of 2 item tuples
        :param ortho: used for linearise and defines which 
                    axis to align along (e.g. align along x with ortho=[1,0,0]), 
                    defaults to [1,1,1]
        :type ortho: list, optional
        :param linearise: make the polymer linear. Good for debugging and
                    visually confirming structure is correct, defaults to False
        :type linearise: bool, optional
        :param beta: beta is used to group monomers into 
                    categories, so that each branch can be separated. Beta is
                    analogous to segments. The tempfactors attribute of all of
                    the atoms in the incoming monomer are set to the value of 
                    beta, defaults to 0
        :type beta: int, optional

        For example: 
            {'P':'CMA','Q':'CA','R':'CN','S':'C'}, with n=1 and 
            nn=2 and joins=[('C','CA')] would produce a polymer of 
            '...-CA1-C1-CA2-C2-...'

        Further notes:
            * n and nn enable branching by specifying from and to connections between monomers and beta is analogous to segments 
            * extend a polymer 'u' by a monomer 'u_', by fitting the backbone atoms (P2, Q2) from the new monomer to (P1,Q1) from the existing residue n, 
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
        
        Note: 
            this function preserves all dummy atoms. Removing dummy atoms 
            during extension causes substantial problems with indexing. They 
            must only be removed once the polymer is fully built. Save your 
            polymer to a .gro or .pdb file with the PDB.save() method, which will strip
            any dummy atoms.

        """

        Q = self.polymer.select_atoms(f"resid {n} and name {names['Q']}").positions[-1]
        S = self.polymer.select_atoms(f"resid {n} and name {names['S']}").positions[-1]

        u_ = monomer # monomer = mdict['path'][monomer]
        u_.atoms.tempfactors=float(beta) 

        u_.residues.resids = nn

        P = u_.select_atoms(f'resid {nn} and name {names["P"]}').positions[0]

        # first, translate u_ by vector P2->P1
        PQ = Q - P

        u_.atoms.translate(PQ)

        # next, rotate around cross product of backbone vectors to align C_n to CN_n+1
        P = u_.select_atoms(f'resid {nn} and name {names["P"]}').positions[0]
        R = u_.select_atoms(f'resid {nn} and name {names["R"]}').positions[0]

        PR = R - P
        PR_n = np.linalg.norm(PR)

        QS =  S - Q
        QS_n = np.linalg.norm(QS)

        theta = degrees(np.arccos(np.dot(PR,QS)/(PR_n * QS_n))) # TODO there are edge cases where this can fail, I think if v1 and v2 are exactly antiparallel. I need to add a checker to resolve that.

        k = np.cross(PR,QS)
        u_r1 = u_.atoms.rotateby(theta,axis=k,point=P)

        if linearise:
            V1= u_.select_atoms('resid '+str(nn)+' and name '+names['V1']).positions[0]
            V2= u_.select_atoms('resid '+str(nn)+' and name '+names['V2']).positions[0]
            V=V2-V1
            V_n=np.linalg.norm(V)
            ortho_n=np.linalg.norm(ortho)
            theta = degrees(np.arccos(np.dot(V,ortho)/(V_n * ortho_n)))

            k=np.cross(V,ortho)
            
            u_r2 = u_r1.atoms.rotateby(theta,axis=k,point=V1)

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
    
    def _split_pol(self,J,J_resid,K,K_resid):
        """
        Given a pair of bonded atoms, uses a graph representation to identify groups of connected atoms on either side of the bond returns atomgroups corresponding to all atoms on either side of the bond.
        
        Use with caution if your polymer contains rings or closed loops, the atoms on either side of the bond must not be connected through other parts of the molecule

        Args:
            J (atom name): the name of the first atom in the bond
            J_resid:       the resid of the first atom in the bond
            K (atom name): the name of the second atom in the bond
            K_resid:       the resid of the second atom in the bond

        Returns two atomgroups, fore and aft
            fore is all atoms connected to J
            aft is all atoms connected to K
        """

        pair = self.polymer.select_atoms(f'(resid {J_resid} and name {J}) or (resid {K_resid} and name {K})' ) 
        bond = self.polymer.atoms.bonds.atomgroup_intersection(pair,strict=True)[0]
        g = nx.Graph()
        g.add_edges_from(self.polymer.atoms.bonds.to_indices()) 
        g.remove_edge(*bond.indices)
        a, _ = (nx.subgraph(g,c) for c in nx.connected_components(g))
        fore=self.polymer.atoms[a.nodes()]
        aft=self.polymer.atoms ^ fore
        return(fore,aft)


    def rotate(self,J,J_resid,K,K_resid,mult=3,step=1):
        """
        Given a pair of bonded atoms J and K within some torsion, uses _split_pol() to identify all atoms connected to J, then rotates them around vector JK by (step * int(360/mult)) degrees, rotating the dihedral centered over J-K by one step.  

        :param J: the name of the first atom in the torsion
        :type J: str
        :param J_resid: the resid of the first atom in the torsion
        :type J_resid: int
        :param K: the name of the second atom in the torsion
        :type K: str
        :param K_resid: the resid of the second atom in the torsion
        :type K_resid: int
        :param mult: the multiplicity of the dihedral centered over J-K, 
                defaults to 3
        :type mult: int, optional
        :param step: how many steps to rotate around the torsion, where one step = 320/mult degrees, 
                defaults to 3
        :type step: int, optional
        """

        fore,_=self._split_pol(J,J_resid,K,K_resid)
        pair = self.polymer.select_atoms(f'(resid {J_resid} and name {J}) or (resid {K_resid} and name {K})' ) 
        bond = self.polymer.atoms.bonds.atomgroup_intersection(pair,strict=True)[0]
        v = bond[1].position - bond[0].position
        o = (fore & bond.atoms)[0]
        rot = step * int(360/mult) # rotate by $step multiplicity units
        fore.rotateby(rot, v, point=o.position)


    def dist(self,J,J_resid,K,K_resid,dummies='X*',backwards_only=True):
        """
        Given a pair of bonded atoms J and K, get minimum distance between
        atoms on one side of a bond, and atoms on the other side of the bond.
        This is useful for detecting overlapping atoms.

        :param J: the name of the first atom in the bond
        :type J: str
        :param J_resid: the resid of the first atom in the bond
        :type J_resid: int
        :param K: the name of the second atom in the bond
        :type K: str
        :param K_resid: the resid of the second atom in the bond
        :type K_resid: int
        :param dummies: the names of dummy atoms, to be excluded from the
                distance calculation, defaults to 'X*'
        :type dummies: str, optional
        :param backwards_only: only consider atoms from residues 1 to
                max([J_resid,K_resid]), and do not consider atoms further
                along the chain. This is useful for solving dihedrals
                algorithmically, as only clashes in the solved region of the
                polymer are considered, defaults to True
        :type backwards_only: bool, optional

        :return: MDAnalysis.analysis.distance_array of the minimum distances
                between atoms on both halves of the bond
        :rtype: mda.analysis.distance_array
        """
        fore,aft=self._split_pol(J,J_resid,K,K_resid)
        trim=''
        if backwards_only:
            maxres=max([J_resid,K_resid])
            trim=f' and (resid 0 to  {maxres})'
        fore_trim=fore.select_atoms(f'(not name {dummies}) {trim}')
        aft_trim=aft.select_atoms(f'(not name {dummies}) {trim}')
        return(distances.distance_array(fore_trim.atoms.positions, aft_trim.atoms.positions).min()) # this is the clash detection

    def shuffle(self,J,J_resid,K,K_resid,dummies='X*',mult=3,cutoff=0.5,clashcheck=False,backwards_only=False):
        """
        Given a pair of bonded atoms J and K, randomly rotates the dihedral
        between them a random amount with clashcheck=True, will check if the
        resulting structure has overlapping atoms, and will undo the rotation
        if so.

        :param J: the name of the first atom in the torsion
        :type J: str
        :param J_resid: the resid of the first atom in the torsion
        :type J_resid: int
        :param K: the name of the second atom in the torsion
        :type K: str
        :param K_resid: the resid of the second atom in the torsion
        :type K_resid: int
        :param dummies: the names of dummy atoms, to be excluded from the
                distance calculation. Passed to :func:`dist`, defaults to 'X*'
        :type dummies: str, optional
        :param mult: the multiplicity of the dihedral centered over J-K, 
                defaults to 3
        :type mult: int, optional
        :param cutoff: maximum interatomic distance for atoms to be considered
                overlapping, defaults to 0.5
        :type cutoff: float, optional
        :param clashcheck: check if the shuffle has introduced a clash, and if
                so reverses the shuffle, defaults to False
        :type clashcheck: bool, optional
        :param backwards_only: passed to :func:`dist` to decide if clash
                checking is done from residue 1 up to max([J_resid,K_resid]) (if True),
                or along the entire polymer (if False)
                defaults to False
        :type backwards_only: bool, optional

        :return: True if clash detected, or False if no clash detected
        :rtype: bool
        """
        step = random.randrange(1,mult) # rotate fore by a random multiplicity
        self.rotate(J,J_resid,K,K_resid,mult,step)
        clash= (self.dist(J,J_resid,K,K_resid,dummies,backwards_only=backwards_only) <= cutoff )
        if clash and clashcheck:
            self.rotate(J,J_resid,K,K_resid,mult,-1 * step,)
        return(clash)

    def dihedral_solver(self,pairlist,dummies='X*',cutoff=0.7,backwards_only=True):
        """
        Converts the shuffled conformation (i.e. after :func:`shuffle`) into
        one without overlapping atoms by resolving dihedrals

        :param pairlist: list of dicts of atom pairs, created with
                :func:`gen_pairlist`, e.g. [{J,J_resid,K,K_resid,mult}]
        :type pairlist: list of dicts of atom pairs
        :param dummies: the names of dummy atoms, to be discluded from the
                distance calculation. Passed to :func:`dist`, defaults to 'X*'
        :type dummies: str, optional
        :param cutoff: maximum interatomic distance for atoms to be considered
                overlapping, defaults to 0.7
        :type cutoff: float, optional
        :param backwards_only: passed to :func:`dist` to decide if clash
                checking is done from residue 1 up to max([J_resid,K_resid]) for the current dihedral (if True),
                or along the entire polymer (if False)
                defaults to True
        :type backwards_only: bool, optional

        :return: True if unable to resolve dihedrals or False if dihedrals all
                resolved and no clashes detected
        :rtype: bool
        """
        steps=len(pairlist)
        tries={x:0 for x in range(0,steps)} # how many steps around the dihedral have we tried? resets to zero if you step backwards
        fails={x:0 for x in range(0,steps)} # how many times have we had to step backwards at this monomer?  the more times, the further back we step 
        # TODO stepback isn't as exhaustive as I'd like it to be.   
        i=0
        failed=False
        done=False
        retry=False
        repeats=0
        while not (done) :
            with tqdm(total=steps) as pbar:
                while i < steps and i >=0:
                    dh=pairlist[i]
                    check=self.dist(J=dh['J'],J_resid=dh['J_resid'],K=dh['K'],K_resid=dh['K_resid'],dummies=dummies,backwards_only=backwards_only)
                    if check > cutoff and not retry:
                        i+=1
                        pbar.update(1)
                    else:
                        retry=False
                        #print(i,'tries',tries[i],'multiplicity:',dh['mult'])
                        if tries[i] >= dh['mult']: # have you tried all steps around the dihedral?
                            if i==0: # yes, and this is the first monomer
                                failed=True
                                done=True
                                i=-1 # force exit from while loops upon failure
                            else: # yes, and this is not the first monomer
                                #print(i,tries[i])
                                retry=True
                                tries[i] = 0 # reset tries for this monomer
                                fails[i]+=1 # monomer has failed
                                pbar.update(-(min(i,fails[i])))
                                i= i-fails[i] # step backwards by number of failures
                                repeats += 1
                        else: # no, there are more tries
                            tries[i] += 1
                            self.rotate(J=dh['J'],J_resid=dh['J_resid'],K=dh['K'],K_resid=dh['K_resid'],mult=dh['mult'])
                done=True
        if failed or i<0: # hard coded to detect failure if you stop at i<=0 because detecting this automatically wasn't working
            print('Could not reach a valid conformation')
            print('Perhaps you should try building a pseudolinear geometry with .extend(linearise=True) or randomising a subset of the dihedrals with shuffler(), and then try solving a conformation again')
            return True
        else:
            return False

    def shuffler(self,pairlist,dummies='X*',cutoff=0.5,clashcheck=False):
        """
        Shuffles each pair of bonded atoms in the provided pairlist with :func:`shuffle`

        :param pairlist: list of dicts of atom pairs, created with
                :func:`gen_pairlist`, e.g. [{J,J_resid,K,K_resid,mult}]
        :type pairlist: list of dicts of atom pairs
        :param dummies: _the names of dummy atoms, to be discluded from the
                distance calculation. Passed to :func:`shuffle`, defaults to 'X*'
        :type dummies: str, optional
        :param cutoff: maximum interatomic distance for atoms to be considered
                overlapping, defaults to 0.5
        :type cutoff: float, optional
        :param clashcheck: check if the shuffle has introduced a clash, and if
                so reverses the shuffle, defaults to False
        :type clashcheck: bool, optional
        """ 
        for dh in tqdm(pairlist):
            self.shuffle(J=dh['J'],J_resid=dh['J_resid'],K=dh['K'],K_resid=dh['K_resid'],mult=dh['mult'],dummies=dummies,clashcheck=clashcheck,cutoff=cutoff)

    def gen_pairlist(self,J,K,first_resid=1,last_resid=999,resid_step=1,same_res=True,mult=3):
        """
        Generates a list of dicts of atom pairs for use in :func:`shuffler` or
        :func:`dihedral_solver`.

        :param J: the name of the first atom in the bond
        :type J: str
        :param K: the name of the second atom in the bond
        :type K: str
        :param first_resid: resid of first residue to include in pairlist,
                defaults to 1
        :type first_resid: int, optional
        :param last_resid: resid of last residue to include in pairlist,
                defaults to 999
        :type last_resid: int, optional
        :param resid_step: interval between each resid, useful for only
                generating a dict in the pairlist for every nth residue,
                defaults to 1
        :type resid_step: int, optional
        :param same_res: if True, paired atoms are bonded within the same
                residue, if False assumes that paired atoms J and K form a
                bond between residue n and residue n+1, defaults to True
        :type same_res: bool, optional
        :param mult: the multiplicity of the dihedral centered over J-K,
                defaults to 3
        :type mult: int, optional

        :return: a list of dicts of atom pairs
        :rtype: list of dicts
        """
        pairlist=[{'J':J,'J_resid':i,'K':K,'K_resid':i+int(not same_res),'mult':mult} for i in range(first_resid,last_resid+1,resid_step)]
        # crude, doesn't actually check if the pairs exist
        # want to extend this to cover more than one dihedral
        return(pairlist)

    def gencomp(mdict, length, fill, middle, count = True, frac = False):
        """
        Deprecated and unsupported, use with caution. 
        Generate the linear conformation of a polymer from a given length, 
        specified monomer composition and a dictionary of monomers. Only used by
        'polyconf_automatic' and should not be used otherwise.
        
        :param mdict: dictionary of resnames making up the 'middle' of the
                polymer (i.e. not terminal residues) with how many of each type
        :type mdict: dict
        :param length: desired length of the polymer (i.e. number of monomers)
        :type length: int
        :param fill: list of monomers which can be used to build the 
                    polymer middle (i.e. are not terminal monomers) if there 
                    are not enough in 'middle' to complete the monomer
        :type fill: list
        :param middle: list of monomers to build the polymer from
        :type middle: list
        :param count: generate absolute composition with counts (e.g. 4 of
                    monomer A and 12 of monomer B), defaults to True
        :type count: bool, optional
        :param frac: generate relative composition with fractions (e.g. 25%
                    monomer A and 75% monomer B), defaults to False
        :type frac: bool, optional

        :return: a list with the sequence of randomly ordered monomers of
                the desired length and composition to build the polymer
        :rtype: list
        """
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

        polyList[0] = polyList[0].split("_")[0] + 'I' # convert first monomer from middle to initiator -> add extra bit onto middle mono
        # TODO: REWORK THIS TO BUILD IN SAME WAY AS POLYTOP!

        polyList[-1] = polyList[-1].split("_")[0] + 'T' # convert last monomer from middle to terminator -> add extra bit onto middle mono
        return(polyList) # list of length items with all monomers in order -> used to build w PDBs
