#!/usr/bin/env python
from .monomer import Monomer
from .PDB import PDB
 
import numpy as np
import pandas as pd
from tqdm import tqdm  # progress bar

import MDAnalysis as mda
from MDAnalysis.analysis import distances
from math import degrees

import random # we'll be using this to randomise monomer positions

import networkx as nx  # we'll be using this to define halves of the molecule during dihedral shuffling

class Polymer:
    """
    Docstrings go brr
    """
    def __init__(self, firstMonomer) -> None:
        # build the first monomer of a polymer
        firstMonomer.residues.resids = 1
        self.first = firstMonomer # is a MDA universe
        self.polymer = self.first
        self.atoms = self.polymer.atoms

    def select_atoms(self, selection):
        return self.polymer.select_atoms(selection)

    def extend(self,monomer,n,nn,names={'P1':'CA','Q1':'C','P2':'CMA','Q2':'CN',},rot=0,joins=[('C','CA')],ortho=[1,1,1],linearise=False,rotate=False,beta=0): # ortho for linearise, e.g. align along x with ortho=[1,0,0]
        # n and nn enable branching by specifying from and to connections between monomers and beta is analogous to segments
        # extend a polymer u by a monomer u_, by fitting the backbone atoms (P2, Q2) from the new monomer to (P1,Q1) from the existing residue n, 
        # then joins the monomer nn to the existing universe with a bond between each pair of atoms in joins
        # ATOM NOMENCLATURE:
        #   P1 and Q1 are atoms in monomer n
        #   P2 and Q2 are dummy atoms in monomer nn, which correspond to the atoms P1 and Q1
        #   R and S are atoms in monomer nn that describe a vector to be aligned to the ortho vector during the linearise fit; 
        #       intended for making straight orthoganal branches for ease of visibility 
        #       this must be chosen carefully so the rotation doe not break stereochemistry or alignment
        # JOINS NOMENCLATURE:
        #   joins contains pairs of atoms to be linked, of the form (X_n,Y_nn)
        #   X_n   is some atom in residue n, the final residue before extension
        #   Y_n+1 is some atom in residue nn, the new residue
        #   I have not tested this script with rings, and I am not sure how it will handle them
        # NB:  this function preserves all dummy atoms. removing dummy atoms during extension causes substantial problems with indexing; you need to remove them later
        # beta is used to group monomers into categories, so that each branch is 

        print("AHHHHHH")
        print(f'resid {n} and name {names["P1"]}')
        print(self.polymer.select_atoms(f'resid {n} and name {names["P1"]}'))
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

        if rotate: # don't worry about it -> use NaK version (prevents atom overlap in adjacent monomers
                R= u_.select_atoms('resid '+str(nn)+' and name '+names['R']).positions[0]
                S= u_.select_atoms('resid '+str(nn)+' and name '+names['S']).positions[0]
                RS=S-R
                RS_n=np.linalg.norm(RS)
                u_r2 = u_r1.atoms.rotateby(rot,axis=RS,point=R)
                u_r1=u_r2

        # readded option to rotate so I can mess with linear chain layouts
        #u_r2 = u_r1.atoms.rotateby(rot,axis=k,point=P1)

        # combine extended polymer into new universe

        new = mda.Merge(self.polymer.atoms, u_r1.atoms)

        # add new bonds linking pairs of atoms (X_n,Y_nn) 

        for pair in joins:
            X = new.select_atoms("resid "+str(n)+" and name "+pair[0]).indices[0]
            Y = new.select_atoms("resid "+str(nn)+" and name "+pair[1]).indices[0]
            new.add_bonds([(X,Y)])

        new.dimensions = list(new.atoms.positions.max(axis=0) + [0.5,0.5,0.5]) + [90]*3
        self.polymer = new.copy()
    
    #REPLACE WITH NAC GENCONF VERSION
    def genconf(self, atomPairs, length, runs = 5, attempts = 20, cutoff = 0.8, fname = 'polymer_conf'):
        # atomPairs is a list of tuples!
        # for bond in self.polymer.bonds: 
            # 
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
                        conf, a, b = self._shuffle(conf,a1, i,a2, i)
                    elif (sign == '+'): # [CA C+]
                        if (i != length):
                            conf, a, b =self._shuffle(conf,a1, i,a2, i+1)
                    else: # [CA C-]
                        if (i != 1):
                            conf, a, b =self._shuffle(conf,a1, i,a2, i-1)

                    #convert a and b from atom indicies to atomgroups
                    aAtoms = conf.atoms[a.nodes()]
                    bAtoms = conf.atoms[b.nodes()]
                    dist=distances.distance_array(aAtoms.atoms.positions, bAtoms.atoms.positions).min() # this is the clash detection
                    if dist >= cutoff:
                        print(f'genconf {i} success')
                        conf.atoms.write(f'{fname}_{i}.pdb')  
                        j+=1
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


    # def genconf(self, # universe containing raw polymer conf with bond information
    # n=5, # number of conformations to generate 
    # cutoff=0.8, # minimum atom/atom distance for clash checker.  This is crude but should be good enough for EM
    # limit=20, # how many times to generate a monomer conformation without clashes, before giving up
    # verbose=False,
    # fname='polymer_conf',
    # ):

    #     errors = {}
    #     e = 0
    #     for rep in range(1,n+1):
    #         conf = self.polymer.copy()  # you need to leave the dummy atoms in, or the bond index doesn't match the atom index
            
    #         # print()
    #         #print(f'\nGenerating conformation {rep} of {n}')
    #         first = True

    #         for i in tqdm(conf.residues.resids,desc=f'Shuffling dihedrals, conformation {rep} of {n}'):
    #             prev = i-1
                
    #             #print(f'\nresid {i}')
    #             #print(conf.select_atoms(f'resid {i}').residues.resnames)

    #             # if not first:
    #             #     # shuffle -CA C bond
    #             #     #print('# shuffle C -CA bond')
    #             #     conf = shuffle(conf,sel=f"(resid {i} and name CA) or (resid {prev} and name C)")
    #             # # shuffle CA C bond
    #             # #print('# shuffle CA C bond')
    #             #done = False
    #             #while not done
    #             conf,clash = self._shuffle(conf, sel=f"resid {i} and name CA C")
    #                 #print(clash)
    #         #        done = clash >= cutoff
    #             # shuffle CA CB bond
    #             #print('# shuffle CA CB bond')
    #             #print(conf.atoms.indices)
    #             #print(len(conf.atoms.indices))

    #             # oh this is failing because of the missing dummy atoms
    #             # ok leave them in for now, cut when you save it
    #             mult=3
    #             done = False
    #             tries = 0
    #             name = conf.select_atoms(f'resid {i}').residues.resnames[0]
    #             while not done:
    #                 conf,clash1 = self._shuffle(conf, sel=f"resid {i} and name CA CB",   mult=mult*2 ) 
    #                 conf,clash2 = self._shuffle(conf, sel=f"resid {i} and name OG1 C1A", mult=mult ) 
    #                 conf,clash3 = self._shuffle(conf, sel=f"resid {i} and name C1A C1B", mult=mult ) 
    #                 clash = min([clash1,clash2,clash3]) # this is quick and dirty, and definitely error prone
    #                 done = clash >= cutoff
    #                 tries += 1

    #                 # if failed, try shuffling entire sidechain

    #                 if tries >= 5 and name in ['P4DM','P4LM','P5DM','P5LM','P4DT','P4LT','P5DT','P5LT','P4DI','P4LI','P5DI','P5LI']:
    #                     if tries == 5 and verbose : print(f'\n!!! ALERT !!!\n\n Clashes remain when shuffling residue {i} {name} after 5 tries.\n Shuffling additional sidechain dihedrals\n')
    #                     conf,clash4 = self._shuffle(conf, sel=f"resid {i} and name C2A C2B", mult=mult*2 ) # bigger space to reduce chance of clashes
    #                     conf,clash5 = self._shuffle(conf, sel=f"resid {i} and name C3A C3B", mult=mult*2 ) # bigger space to reduce chance of clashes
    #                     conf,clash6 = self._shuffle(conf, sel=f"resid {i} and name C4A C4B", mult=mult*2 ) # bigger space to reduce chance of clashes
    #                     clash = min([clash,clash4,clash5,clash6])
    #                     done = clash >= cutoff

    #                 if tries >= 10 and not first:
    #                     conf,clash7 = self._shuffle(conf,sel=f"(resid {prev} and name CA) or (resid {prev} and name C)", mult=mult*2) # try shuffling the -C CA bond too
    #                     conf,clash8 = self._shuffle(conf,sel=f"(resid {prev} and name C) or (resid {i} and name CA)", mult=mult*2) # try shuffling the -C CA bond too
    #                     # generally only do this to resolve clashes; I've found that shuffling both -C to CA and CA to C in every monomer leads to tightly packed structures with many clashes

    #                     clash = min([clash,clash7,clash8])
    #                     done = clash >= cutoff

    #                 if tries == 20:
    #                     mult=6 # double dihedral multiplicity to increase search space for conformations without clashes

    #                 if tries >= limit: 
    #                     print(f'\n\n!!!!!!!!!!!\n! WARNING !\n!!!!!!!!!!!\n\n Clashes remain when shuffling residue {i} {name} after {limit} tries.\n Atoms remain within {cutoff} angstroms.\n Continuing with this geometry.\n You may wish to check for overlapping atoms.\n')
    #                     errors[e]={'conf':rep,'resid':i,'name':name}
    #                     e+=1
    #                     done = True

    #             # if name in ['P5DM','P5LM','P5DT','P5LT','P5DI','P5LI']:
    #             #     conf,clash = shuffle(conf, sel=f"resid {i} and name C5A C5B", mult=6 ) # bigger space to reduce chance of clashes

    #                 #print(clash)
    #             #    done = clash >= cutoff

    #             first=False # used for if you want to do -CA C shuffling TO DO: What is this???

    #         confSaver = PDB(conf)
    #         confSaver.cleanup()
    #         confSaver.save(fname = f'{fname}_{rep}.gro')

    #     print('\n\nFinished generating conformations\n\n')

    #     if e > 0:
    #         print('Clashes log:\n')
    #         for v in errors.values():
    #             conf,resid,name = v.values()
    #             print(f'When generating confromation {conf}, clash at residue {resid} {name} ')

   
    def _shuffle(self, u, a1, a1resid, a2, a2resid, mult=3):
        # based on a tutorial by richard j gowers; http://www.richardjgowers.com/2017/08/14/rotating.html
        pair = u.select_atoms(f'(resid {a1resid} and name {a1}) or (resid {a2resid} and name {a2})' ) 
        #print(pair.atoms)
        bond = u.atoms.bonds.atomgroup_intersection(pair,strict=True)[0]
        #print(bond.atoms)
        g = nx.Graph()
        g.add_edges_from(u.atoms.bonds.to_indices()) # get entire residue as a graph
        g.remove_edge(*bond.indices)     # remove the bond from the graph
        # unpack the two unconnected graphs
        a, b = (nx.subgraph(g,c) for c in nx.connected_components(g))
        # call the graph without the CA atom the 'head', this will be rotated
        fore_nodes = a
        fore = u.atoms[fore_nodes.nodes()]
        v = bond[1].position - bond[0].position
        o = (fore & bond.atoms)[0]
        rot = random.randrange(0,mult)*int(360/mult) # rotate fore by a random multiplicity
        fore.rotateby(rot, v, point=o.position)
        return(u, a, b)
    
    # def _shuffle(self, polymerCopy,sel='name CA C',mult=3): # doesn't modify self polymer incase genconf cannot resolve conformation
    #     # based on a tutorial by richard j gowers; http://www.richardjgowers.com/2017/08/14/rotating.html
    #     pair = polymerCopy.select_atoms(sel)
    #     resmin = pair.residues.resids.min()
    #     resmax = pair.residues.resids.max()
    #     #print(pair.atoms)
    #     print(polymerCopy.select_atoms('name CIA'))
    #     CIA = polymerCopy.select_atoms('name CIA')[0].index # use the initiator alpha carbon to define the fore group
    #     bond = polymerCopy.atoms.bonds.atomgroup_intersection(pair,strict=True)[0]
    #     #print(bond.atoms)
    #     g = nx.Graph()
    #     g.add_edges_from(polymerCopy.atoms.bonds.to_indices()) # get entire polymer as a graph
    #     g.remove_edge(*bond.indices)     # remove the bond from the graph
    #     i=0
    #     for c in nx.connected_components(g):    
    #         #print(i,c)
    #         i+=1
    #     # unpack the two unconnected graphs
    #     a, b = (nx.subgraph(g,c) for c in nx.connected_components(g))
    #     # call the graph with the initiator CA atom the 'head'
    #     fore_nodes = a if CIA in a.nodes() else b
    #     fore = polymerCopy.atoms[fore_nodes.nodes()]
    #     aft = polymerCopy.atoms ^ fore
    #     v = bond[1].position - bond[0].position
    #     o = (aft & bond.atoms)[0]
    #     rot = random.randrange(0,mult)*int(360/mult) # rotate by a random multiplicity
    #     aft.rotateby(rot, v, point=o.position)
    #     forechk = fore.select_atoms(f'(resid 1 to {resmin}) and (not name CP CQ CN CMA)') # dummy atoms can't clash
    #     aftchk = aft.select_atoms(f'(resid 1 to {resmax}) and (not name CP CQ CN CMA)') # dummy atoms can't clash
    #     dist = round(distances.distance_array(forechk.atoms.positions, aftchk.atoms.positions).min(),1)  # minimum distance between fore and aft atoms, not used for now
    #     #print(dist)
    #     return(polymerCopy, dist)
    
    def gencomp(mdict, length, fill, middle, count = True, frac = False): # args are dictionary of resnames in the middle with how many of each type, the total monomer count and a list w random monomers to "bulk" out with (same as from mdict)
        polycomp = []

        if count: # use a count - i.e. 75 of thing A
            for m in middle:
                rcount = mdict['count'][m]
                for i in range(rcount):
                    polycomp += [m] # build a list with one of each monomer unit present in the final polymer

        elif frac: # use a fraction - i.e. 50% of thing A
            for m in middle:
                rfrac = int(length * mdict['frac'][m]) # fraction of total length that is monomer x
                for i in range(rfrac):
                    polycomp += [m] # build a list with one of each monomer unit present in the final polymer
        # then randomise positions of monomers'

        # 'polycomp' is ordered list of X number of monomers -> will add all of monomer A required, before doing all of B, etc.

        filler = [x for x in fill]
        # print(filler)
        # print(fill)
        while len(polycomp) < length:
            sel = random.sample(filler,1)
            polycomp += sel # if the polymer is shorter than the desired length, continue to fill from the given list 
            filler = [x for x in filler if not x in sel]
            if len(filler) == 0: 
                filler = [x for x in fill]

        polyList=random.sample(polycomp,length) # reorder the list of monomers randomly

        polyList[0] = polyList[0][:-1] + 'I' # convert first monomer from middle to initiator -> add extra bit onto middle mono
        # REWORK THIS TO BUILD IN SAME WAY AS POLYTOP!

        polyList[-1] = polyList[-1][:-1] + 'T' # convert last monomer from middle to terminator -> add extra bit onto middle mono
        return(polyList) # list of length items with all monomers in order -> used to build w PDBs **FIND WAY TO YEET INTO POLYTOP!!**
