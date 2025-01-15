#!/usr/bin/env python 
from .Monomer import Monomer
from .Polymer import Polymer
from .PDB import PDB
from polytop.polytop_automatic import Automatic

import numpy as np
import pandas as pd
from tqdm import tqdm  # progress bar

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--name", default='polymer', help='system name string, used for filenames')
parser.add_argument("--nconfs", default=3, )
parser.add_argument('--count', action='store_true', help='explicit count of how many of each type of monomer to include')
parser.add_argument('--frac', action='store_true', help='fraction of polymer made up of each type of monomer')
parser.add_argument('--length', type=int,required=True, help='total polymer length')
parser.add_argument("--monomers",type=str,required=True,default='monomers.csv',help='path to csv describing monomers, expected columns are "resname", "pdb path", "itp path", "position", "fill", and at least one of "count" or "frac"') 
parser.add_argument('--shuffles', type=int,default=20, help='max number of shuffles before restarting conformation generation')
parser.add_argument('--rotate', metavar='N', type=str, nargs='*', required=True,
                    help='which bond pairs to shuffle, e.g. "CA C CA CB" to shuffle CA-C and CA-CB')
parser.add_argument('--joiners', metavar='N', type=str, nargs='*', required=True,
                    help='which atoms will be the corresponding dummies for polyconf extend and the corresponding "residue atoms" for the "monomer atoms" specified in --junctions, e.g. "CMA CN", where the first junction atom must be bonded to the first dummy atom, and the second junction atom to the second dummy atom')
parser.add_argument('--junctions', metavar='N', type=str, nargs='*', required=True,
                    help='which atoms will be "monomer atoms" for junctions of monomers joined together, define in order as "from" atom, then "to" atom, e.g. "C CA" makes bond from C_i to CA_i+1')

# assumes monomer resnames are four letter codes, with the final letter either M for middle, T for terminator, or I for initiator)
# count = 
# frac = fraction of polymer that is made up of these monomers
# path = path to coordinate file with bond information, eg pdb file with connect records
# position = initial, middle, terminal
# fill = after polymer is extended based on count or frac, gencomp checks if the desired  length is required
#           if the polymer has not reached the desired length, extends by choosing randomly from monomers where fill = True

def main():
    print("\nWARNING")
    print("This automatic linear, united-atom polymer builder is newly developed and may not "
          "construct polymer topologies and coordinates as expected, please "
          "manually check all output for correctness before submitting for MD\n\n")
    args = parser.parse_args()

    fname = '_'.join(args.name.split(' '))

    df = pd.read_csv(args.monomers) # has paths to all ITP and PDB files
    df.set_index('resname',inplace=True)
    if args.count:
        df['count'] = df['count'].astype(int) 
    elif args.frac:
        df['frac'] = df['frac'].astype(float)
    else:
        raise(ValueError, "Must specify either '--count' or '--frac' for polymer composition")
    df['fill'] = df['fill'].astype(int).astype(bool)
    mdict = df.to_dict()

    # build iterator lists 

    monomers = [x for x in mdict['pdb path'].keys()]

    middle = [x for x in monomers if mdict['position'][x] == 'middle']    
    initial = [x for x in monomers if mdict['position'][x] == 'initial']
    terminal = [x for x in monomers if mdict['position'][x] == 'terminal']

    middle.sort()
    initial.sort()
    terminal.sort()

    monomers = [*initial,*middle,*terminal]

    fill = [x for x in monomers if mdict['fill'][x]]

    polyList = Polymer.gencomp(mdict, args.length, fill, middle, count = args.count, frac = args.frac)
    print("Generated polymer composition monomer order:")
    print(polyList)

    print('\nPolymer composition is:')
    for l in [initial, middle, terminal]:
        for m in l:
            print(m,':',len([x for x in polyList if x == m]))

    print(f"\nPolymer building by joining atom {args.junctions[0]} to atom {args.junctions[1]}")

    # has paths to all ITP and PDB files
    itp_monomers = [x for x in mdict['itp path'].keys()]

    polymerITP = Automatic(polyList, mdict['itp path'], args.length, itp_monomers, args.junctions, args.joiners)
    print("\n\nBuilding polymer topology with PolyTop...\n")
    polymerITP.build(outputName=args.name)
    print(f"Saved polymer topology as '{fname}.itp'\n")



    # Generate and save linear polymer
    j_a = args.junctions
    d_a = args.joiners
    polymer = Polymer(Monomer(mdict['pdb path'][polyList[0]]))
    for i in tqdm(range(len(polyList)),desc='Building initial polymer geometry with PolyConf'):
        if i >= 1:
            names = {'P1':d_a[1],'Q1':j_a[0],'P2':j_a[1],'Q2':d_a[0], 'R':args.junctions[1], 'S':args.junctions[0]}
            j_t = (j_a[0], j_a[1])
            joins = [j_t]
            polymer.extend(Monomer(mdict['pdb path'][polyList[i]]),n=i,nn=i+1,names=names, joins=joins, linearise=True)
            # rename the dummy atoms to 'X' to flag for removal
            for j, name in enumerate(args.joiners):
                polymer.renamer(i+j, name, nameout='X')
        
    polymerSaver = PDB(polymer)
    polymerSaver.cleanup()
    polymerSaver.save(fname = f'{fname}_linear', selectionString=None)
    print(f"\n\nSaved linear polymer geometry as '{fname}_linear.pdb'\n\n")


    # Copy linear polymer, generate different conformations without steric clashes and save
    CN=polymer.gen_pairlist(a1=args.rotate[0],a2=args.rotate[1],first_resid=1,last_resid=args.length,mult=3)

    for i in range(args.nconfs):
        failed = True
        polymerToShuffle = polymer.copy()
        for j in range(args.shuffles):
            failed = polymerToShuffle.dihedral_solver(CN,dummy=args.joiners,cutoff=0.7)
            if failed == False:
                print(f"successfully generated conformation {i+1} without clashes")
                Saver = PDB(polymerToShuffle)
                Saver.cleanup() # center in box
                Saver.save(dummyAtoms=args.joiners,fname=f'{args.name}_solved{i+1}')
                print(f"Saved polymer geometry {i+1} as '{fname}_solved{i+1}.pdb'\n\n")
                break
        if failed == True:
            print(f"Unable to generate conformation number {i+1} - moving onto the next conformation\n\n")