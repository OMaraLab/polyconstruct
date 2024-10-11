#!/usr/bin/env python 
from polyconf.monomer import Monomer
from polyconf.polymer import Polymer
from polyconf.PDB import PDB
# from polytop import Automatic

import numpy as np
import pandas as pd
from tqdm import tqdm  # progress bar

import MDAnalysis as mda
from MDAnalysis.analysis import distances
from math import degrees

import random # we'll be using this to randomise monomer positions

import networkx as nx  # we'll be using this to define halves of the molecule during dihedral shuffling
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--name", default='polymer', help='system name string, used for filenames')
parser.add_argument("--nconfs", default=3)
parser.add_argument('--count', action='store_true') # specifies monomers are given by count
parser.add_argument('--frac', action='store_true') # specifies monomers are as fractions, requires a total length
parser.add_argument('--length', type=int,required=True) # total polymer legnth
parser.add_argument("--monomers",type=str,default='monomers.csv',help='path to csv describing monomers, expected columns are \'resname\', \'path\', \'position\', \'fill\', and at least one of \'count\' and \'frac\'') 
parser.add_argument('--shuffles', type=int,default=20) # specifies monomers are given by count

# assumes monomer resnames are four letter codes, with the final letter either M for middle, T for terminator, or I for initiator)
# count = explicit count of how many to include
# frac = fraction of polymer that is made up of these monomers
# path = path to coordinate file with bond information, eg pdb file with connect records
# position = initial, middle, terminal
# fill = after polymer is extended based on count or frac, gencomp checks if the desired  length is required
#           if the polymer has not reached the desired length, extends by choosing randomly from monomers where fill = True

args = parser.parse_args()

fname = '_'.join(args.name.split(' '))
nconfs = int(args.nconfs)

df = pd.read_csv(args.monomers) # has paths to all ITP and PDB files
df.set_index('resname',inplace=True)
if args.count:
    df['count'] = df['count'].astype(int) 
if args.frac:
    df['frac'] = df['frac'].astype(float) 
df['fill'] = df['fill'].astype(int).astype(bool)
mdict = df.to_dict()

# build iterator lists 

monomers = [x for x in mdict['path'].keys()]

middle = [x for x in monomers if mdict['position'][x] == 'middle']    
initial = [x for x in monomers if mdict['position'][x] == 'initial']
terminal = [x for x in monomers if mdict['position'][x] == 'terminal']

middle.sort()
initial.sort()
terminal.sort()

monomers = [*initial,*middle,*terminal]

fill = [x for x in monomers if mdict['fill'][x]]
print(fill)

polyList = Polymer.gencomp(mdict, args.length, fill, middle, count = args.count, frac = args.frac)
print(polyList)

# df2 = pd.read_csv(args.monomers) # has paths to all ITP and PDB files
# ITPLocations = df2['ITPpath'].tolist() #TO DO: create test CSV for this!
# junctions = df2['junctions'].tolist() # get junctions from polyconf?

# polymerITP = Automatic(polyList, args.length, ITPLocations)
# polymerITP.build()
# TO DO: create new class to automate polytop workflow!!
# now, use the composition to build the polymer

print('Polymer composition generated\n')
for l in [initial, middle, terminal]:
    for m in l:
        print(m,':',len([x for x in polyList if x == m]))
    print()

polymer = Polymer(Monomer(mdict['path'][polyList[0]]))
for i in tqdm(range(len(polyList)),desc='Building initial polymer geometry'):
    print(mdict['path'][polyList[i]])
    

    if i >= 1:
        # print(polyList[i])
        # print(mdict['path'][polyList[i]])
        
        
        #print( pol.residues.resids)
        #print(pol.select_atoms('resid 1 and name CA').positions[0])
    # else:
        rot = (i*45)%360 # not used anymore; replaced with genconf shuffling
        polymer.extend(Monomer(mdict['path'][polyList[i]]),n=i,nn=i+1,rot=rot)


polymerSaver = PDB(polymer)
polymerSaver.cleanup()
polymerSaver.save(fname = f'{fname}_linear.gro', selectionString=None)

# polymer.genconf(n=nconfs,fname=fname,verbose=False,limit=args.shuffles)

atomPairs = [('CA', 'C'), ('CA', 'CB')] # specified by user
# for i in range(polymer.bonds):
#     for atom in polymer.atoms:
#         if atom in polymer.bonds[i]:
#             atomPairs[i] = (atom, polymer.bonds[i].partner(atom))

polymer.genconf(atomPairs, length = args.length, runs = 5, cutoff = 0.8, fname = 'polymer_conf')