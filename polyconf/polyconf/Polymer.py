#!/usr/bin/env python
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
    
    def copy(self) -> mda.Universe:
        """
        Use MDAnalysis Merge to create a new MDAnalysis Universe that is an 
        exact copy of this one containing the polymer. Deepcopy is not used as
        it can be problematic with open file sockets.

        :return: A new MDAnalysis Universe that is an exact copy of the polymer
        :rtype: mda.Universe
        """
        new_u = mda.Merge(self.polymer.atoms)
        return Polymer(Monomer.monomer_from_u(new_u), keep_resids=True)
        # return Polymer(Monomer.monomer_from_u(new_u), keep_resids=True)
    
    def renamer(self, resid, namein, nameout='X'):
        """
        Change selected atom names to a new name like X1, X2, X3. Intended to
        flag dummy atoms for removal. Selected atoms are given a basename, 
        e.g. 'X' defined by the nameout argument, as well as a number. 

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
        Generates a resid that is one larger than the highest existing resid in
        the polymer. Excellent for finding what the number of the next residue
        should be set to, for the nn argument to 'extend()'.

        :return: resid one greater than the polymer's current highest resid
        :rtype: int
        """
        nn = max(self.polymer.residues.resids) + 1
        return nn


    def maxresid(self):
        """
        Returns maximum resid in the polymer
        
        :return: the polymer's current highest resid
        :rtype: int
        """
        n = max(self.polymer.residues.resids)
        return n


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
        :param names: four single 'key:value' pairs with keys
                    P1, Q1, P2 and Q2 for the four dummy atoms that are 
                    'overlaid' during polymer extension. P1 is the real atom in
                    the polymer while P2 is its corresponding dummy atom in the
                    incoming monomer. Similarly, Q1 is the real atom in the 
                    polymer while Q2 is the corresponding dummy atom in the 
                    incoming polymer. Such that P1 and Q1 and bonded, and P2 
                    and Q2 have an equivalent bond. The result will be a new 
                    bond between Q1 and the next atom beside Q2 that is NOT P2.
        :type names: dict
        :param joins: defines which two atoms of residue
                    'n' (first value) and residue 'nn' (second value) will 
                    be bonded
        :type joins: list of a 2 item tuple
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
            {'P1':'CA','Q1':'C','P2':'CMA','Q2':'CN'}, with n=1 and 
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
        
        Note: 
            this function preserves all dummy atoms. Removing dummy atoms 
            during extension causes substantial problems with indexing. They 
            must only be removed once the polymer is fully built. Save your 
            polymer to a .gro or .pdb file with the PDB.save() method, which will strip
            any dummy atoms.

        """
        # TODO I renamed the atoms in the manuscript, so I guess I gotta fix that here
        # cool.
        # love that for me.

        P1 = self.polymer.select_atoms(f'resid {n} and name {names["P1"]}').positions[-1]
        Q1 = self.polymer.select_atoms(f'resid {n} and name {names["Q1"]}').positions[-1]

        u_ = monomer # monomer = mdict['path'][monomer]
        u_.atoms.tempfactors=float(beta) 

        u_.residues.resids = nn

        P2 = u_.select_atoms(f'resid {nn} and name {names["P2"]}').positions[0]

        # first, translate u_ by vector P2->P1
        T = P1 - P2

        u_.atoms.translate(T)

        # next, rotate around cross product of backbone vectors to align C_n to CN_n+1
        P2 = u_.select_atoms(f'resid {nn} and name {names["P2"]}').positions[0]
        Q2 = u_.select_atoms(f'resid {nn} and name {names["Q2"]}').positions[0]

        v1 = Q2 - P2
        v1_n = np.linalg.norm(v1)

        v2 =  Q1 - P1
        v2_n = np.linalg.norm(v2)

        theta = degrees(np.arccos(np.dot(v1,v2)/(v1_n * v2_n))) # TODO there are edge cases where this can fail, I think if v1 and v2 are exactly antiparallel. I need to add a checker to resolve that.

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

            # R= u_.select_atoms('resid '+str(nn)+' and name '+names['R']).positions[0]
            # S= u_.select_atoms('resid '+str(nn)+' and name '+names['S']).positions[0]
            # RS=S-R
            # RS_n=np.linalg.norm(RS)
            # ortho_n=np.linalg.norm(ortho)
            # check = degrees(np.arccos(np.dot(RS,ortho)/(RS_n * ortho_n)))
            # if abs(check)>1: print('linear check =' ,check)

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
    
    def _split_pol(self,a1,a1_resid,a2,a2_resid):
        """
        Given a pair of bonded atoms, uses a graph representation to identify groups of connected atoms on either side of the bond returns atomgroups corresponding to all atoms on either side of the bond.
        
        Use with caution if your polymer contains rings or closed loops, the atoms on either side of the bond must not be connected through other parts of the molecule

        Args:
            a1 (atom name): the name of the first atom in the bond
            a1_resid:       the resid of the first atom in the bond
            a2 (atom name): the name of the second atom in the bond
            a2_resid:       the resid of the second atom in the bond

        Returns two atomgroups, fore and aft
            fore is all atoms connected to a1
            aft is all atoms connected to a2
        """

        pair = self.polymer.select_atoms(f'(resid {a1_resid} and name {a1}) or (resid {a2_resid} and name {a2})' ) 
        bond = self.polymer.atoms.bonds.atomgroup_intersection(pair,strict=True)[0]
        g = nx.Graph()
        g.add_edges_from(self.polymer.atoms.bonds.to_indices()) 
        g.remove_edge(*bond.indices)
        a, _ = (nx.subgraph(g,c) for c in nx.connected_components(g))
        fore=self.polymer.atoms[a.nodes()]
        aft=self.polymer.atoms ^ fore
        return(fore,aft)


    def _rotate(self,a1,a1_resid,a2,a2_resid,mult=3,step=1):
        """
        Given a pair of bonded atoms a1 and a2, uses _split_pol() to identify all atoms connected to a1, then rotates them around the bond by (step * int(360/mult)) degrees, chanding the state of the dihedral centered over a1-b1.  

        Args:
            a1 (atom name): the name of the first atom in the bond
            a1_resid:       the resid of the first atom in the bond
            a2 (atom name): the name of the second atom in the bond
            a2_resid:       the resid of the second atom in the bond
            mult:           the multiplicity of the dihedral centered over a1-b1
            step:           how many steps around the dihedral to rotate
        """

        fore,_=self._split_pol(a1,a1_resid,a2,a2_resid)
        pair = self.polymer.select_atoms(f'(resid {a1_resid} and name {a1}) or (resid {a2_resid} and name {a2})' ) 
        bond = self.polymer.atoms.bonds.atomgroup_intersection(pair,strict=True)[0]
        v = bond[1].position - bond[0].position
        o = (fore & bond.atoms)[0]
        rot = step * int(360/mult) # rotate by $step multiplicity units
        fore.rotateby(rot, v, point=o.position)


    def dist(self,a1,a1_resid,a2,a2_resid,dummy='X*',backwards_only=True):
        """
        Given a pair of bonded atoms a1 and a2, get minimum distance between
        atoms on one side of a bond, and atoms on the other side of the bond.
        This is useful for detecting overlapping atoms.

        :param a1: the name of the first atom in the bond
        :type a1: str
        :param a1_resid: the resid of the first atom in the bond
        :type a1_resid: int
        :param a2: the name of the second atom in the bond
        :type a2: str
        :param a2_resid: the resid of the second atom in the bond
        :type a2_resid: int
        :param dummy: the names of dummy atoms, to be discluded from the
                distance calculation, defaults to 'X*'
        :type dummy: str, optional
        :param backwards_only: only consider atoms from residues 1 to
                max([a1_resid,a2_resid]), and do not consider atoms further
                along the chain. This is useful for solving dihedrals
                algorithmically, as only clashes in the solved region of the
                polymer are considered, defaults to True
        :type backwards_only: bool, optional

        :return: MDAnalysis.analysis.distance_array of the minimum distances
                between atoms on both halves of the bond
        :rtype: mda.analysis.distance_array
        """
        fore,aft=self._split_pol(a1,a1_resid,a2,a2_resid)
        trim=''
        if backwards_only:
            maxres=max([a1_resid,a2_resid])
            trim=f' and (resid 0 to  {maxres})'
        fore_trim=fore.select_atoms(f'(not name {dummy}) {trim}')
        aft_trim=aft.select_atoms(f'(not name {dummy}) {trim}')
        return(distances.distance_array(fore_trim.atoms.positions, aft_trim.atoms.positions).min()) # this is the clash detection

    def shuffle(self,a1,a1_resid,a2,a2_resid,dummy='X*',mult=3,cutoff=0.5,clashcheck=False,backwards_only=False):
        """
        Given a pair of bonded atoms a1 and a2, randomly rotates the dihedral
        between them a random amount with clashcheck=True, will check if the
        resulting structure has overlapping atoms, and will undo the rotation
        if so.

        :param a1: the name of the first atom in the bond
        :type a1: str
        :param a1_resid: the resid of the first atom in the bond
        :type a1_resid: int
        :param a2: the name of the second atom in the bond
        :type a2: str
        :param a2_resid: the resid of the second atom in the bond
        :type a2_resid: int
        :param dummy: the names of dummy atoms, to be discluded from the
                distance calculation. Passed to :func:`dist`, defaults to 'X*'
        :type dummy: str, optional
        :param mult: the multiplicity of the dihedral centered over a1-b1, 
                defaults to 3
        :type mult: int, optional
        :param cutoff: maximum interatomic distance for atoms to be considered
                overlapping, defaults to 0.5
        :type cutoff: float, optional
        :param clashcheck: check if the shuffle has introduced a clash, and if
                so reverses the shuffle, defaults to False
        :type clashcheck: bool, optional
        :param backwards_only: passed to :func:`dist` to decide if clash
                checking is done along the entire polymer (if True) or just
                from residue 1 up to max([a1_resid,a2_resid]) (if False),
                defaults to False
        :type backwards_only: bool, optional

        :return: True if clash detected, or False if no clash detected
        :rtype: bool
        """
        step = random.randrange(1,mult) # rotate fore by a random multiplicity
        self._rotate(a1,a1_resid,a2,a2_resid,mult,step)
        clash= (self.dist(a1,a1_resid,a2,a2_resid,dummy,backwards_only=backwards_only) <= cutoff )
        if clash and clashcheck:
            self._rotate(a1,a1_resid,a2,a2_resid,mult,-1 * step,)
        return(clash)

    def dihedral_solver(self,pairlist,dummy='X*',cutoff=0.7):
        """
        Converts the shuffled conformation (i.e. after :func:`shuffle`) into
        one without overlapping atoms by resolving dihedrals

        :param pairlist: list of dicts of atom pairs, created with
                :func:`gen_pairlist`, e.g. [{a1,a1_resid,a2,a2_resid,mult}]
        :type pairlist: list of dicts of atom pairs
        :param dummy: the names of dummy atoms, to be discluded from the
                distance calculation. Passed to :func:`dist`, defaults to 'X*'
        :type dummy: str, optional
        :param cutoff: maximum interatomic distance for atoms to be considered
                overlapping, defaults to 0.7
        :type cutoff: float, optional

        :return: True if unable to resolve dihedrals or False if dihedrals all
                resolved and no clashes detected
        :rtype: bool
        """
        # TODO stepback isn't working right, it's not stepping back far enough
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
                    check=self.dist(a1=dh['a1'],a1_resid=dh['a1_resid'],a2=dh['a2'],a2_resid=dh['a2_resid'],dummy=dummy)
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
                            self._rotate(a1=dh['a1'],a1_resid=dh['a1_resid'],a2=dh['a2'],a2_resid=dh['a2_resid'],mult=dh['mult'])
                done=True
        if failed or i<0: # hard coded to detect failure if you stop at i<=0 because detecting this automatically wasn't working
            print('Could not reach a valid conformation')
            print('Perhaps you should try building a pseudolinear geometry with .extend(linearize=True) or randomising all dihedrals with shuffler(), and then try solving a conformation again') 
        return(failed or i<0)

    def shuffler(self,pairlist,dummy='X*',cutoff=0.5,clashcheck=False):
        """
        Shuffles each pair of bonded atoms in the provided pairlist with :func:`shuffle`

        :param pairlist: list of dicts of atom pairs, created with
                :func:`gen_pairlist`, e.g. [{a1,a1_resid,a2,a2_resid,mult}]
        :type pairlist: list of dicts of atom pairs
        :param dummy: _the names of dummy atoms, to be discluded from the
                distance calculation. Passed to :func:`shuffle`, defaults to 'X*'
        :type dummy: str, optional
        :param cutoff: maximum interatomic distance for atoms to be considered
                overlapping, defaults to 0.5
        :type cutoff: float, optional
        :param clashcheck: check if the shuffle has introduced a clash, and if
                so reverses the shuffle, defaults to False
        :type clashcheck: bool, optional
        """ 
        for dh in tqdm(pairlist):
            self.shuffle(a1=dh['a1'],a1_resid=dh['a1_resid'],a2=dh['a2'],a2_resid=dh['a2_resid'],mult=dh['mult'],dummy=dummy,clashcheck=clashcheck,cutoff=cutoff)

    def gen_pairlist(self,a1,a2,first_resid=1,last_resid=999,resid_step=1,same_res=True,mult=3):
        """
        Generates a list of dicts of atom pairs for use in :func:`shuffler` or
        :func:`dihedral_solver`.

        :param a1: the name of the first atom in the bond
        :type a1: str
        :param a2: the name of the second atom in the bond
        :type a2: str
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
                residue, if False assumes that paired atoms a1 and a2 form a
                bond between residue n and residue n+1, defaults to True
        :type same_res: bool, optional
        :param mult: the multiplicity of the dihedral centered over a1-b1,
                defaults to 3
        :type mult: int, optional

        :return: a list of dicts of atom pairs
        :rtype: list of dicts
        """
        pairlist=[{'a1':a1,'a1_resid':i,'a2':a2,'a2_resid':i+int(not same_res),'mult':mult} for i in range(first_resid,last_resid+1,resid_step)]
        # crude, doesn't actually check if the pairs exist
        # want to extend this to cover more than one dihedral
        return(pairlist)

    def gencomp(mdict, length, fill, middle, count = True, frac = False):
        """
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
