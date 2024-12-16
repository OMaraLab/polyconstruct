#!/usr/bin/env python 
from .Monomer import Monomer
from .Polymer import Polymer
from .PDB import PDB
import polytop
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
                    help='which atoms will be the corresponding dummies for polyconf extend, e.g. "CMA CN"')
parser.add_argument('--junctions', metavar='N', type=str, nargs='*', required=True,
                    help='which atoms will be "monomer atoms" for junctions of monomers joined together, e.g. "C CA"')
parser.add_argument('--dummies', metavar='N', type=str, nargs='*', required=True,
                    help='which atoms will be the corresponding "residue atoms" for the "monomer atoms" specified in --junctions, e.g. "CP CN"')

# assumes monomer resnames are four letter codes, with the final letter either M for middle, T for terminator, or I for initiator)
# count = 
# frac = fraction of polymer that is made up of these monomers
# path = path to coordinate file with bond information, eg pdb file with connect records
# position = initial, middle, terminal
# fill = after polymer is extended based on count or frac, gencomp checks if the desired  length is required
#           if the polymer has not reached the desired length, extends by choosing randomly from monomers where fill = True

def main():
    print("WARNING")
    print("this automatic linear polymer builder is newly developed and may not "
          "construct polymer topologies and coordinates as expected, please "
          "manually check all output for correctness before submitting for MD")
    args = parser.parse_args()

    fname = '_'.join(args.name.split(' '))
    nconfs = int(args.nconfs)

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

    # has paths to all ITP and PDB files
    itp_monomers = [x for x in mdict['itp path'].keys()]

    polymerITP = Automatic(polyList, mdict['itp path'], args.length, itp_monomers, args.junctions, args.dummies)
    print("\nBuilding polymer topology")
    polymerITP.build(outputName=args.name)

    print('\nPolymer composition:')
    for l in [initial, middle, terminal]:
        for m in l:
            print(m,':',len([x for x in polyList if x == m]))

    # Generate and save linear polymer
    j_a = args.junctions
    d_a = args.joiners
    polymer = Polymer(Monomer(mdict['pdb path'][polyList[0]]))
    for i in tqdm(range(len(polyList)),desc='Building initial polymer geometry'):
        if i >= 1:
            rot = (i*45)%360 # not used anymore; replaced with genconf shuffling
            # names = {'P1':'CA','Q1':'C','P2':'CMA','Q2':'CN'}
            names = {'P1':j_a[1],'Q1':j_a[0],'P2':d_a[0],'Q2':d_a[1], 'R':args.junctions[0], 'S':args.junctions[1]}
            j_t = (j_a[0], j_a[1])
            joins = [j_t]
            polymer.extend(Monomer(mdict['pdb path'][polyList[i]]),n=i,nn=i+1,names=names, joins=joins, linearise=True)

    polymerSaver = PDB(polymer)
    polymerSaver.cleanup()
    polymerSaver.save(args.dummies, fname = f'{fname}_linear', selectionString=None)


    # Copy linear polymer, generate different conformations without steric clashes and save
    CN=polymer.gen_pairlist(a1=args.rotate[0],a2=args.rotate[1],first_resid=1,last_resid=args.length,mult=3)

    for i in range(args.nconfs):
        success = True
        polymerToShuffle = polymer.copy()
        for j in range(args.shuffles):
            success = polymerToShuffle.dihedral_solver(CN,dummy=args.joiners,cutoff=0.7)
            if success == False:
                print(f"Successfully generated conformation {i+1} without clashes")
                Saver = PDB(polymerToShuffle)
                Saver.cleanup() # center in box
                Saver.save(dummyAtoms=args.joiners,fname=f'{args.name}_solved{i+1}')
                break
        if success == True:
            print(f"Unable to generate conformation number {i+1} - moving onto the next conformation")