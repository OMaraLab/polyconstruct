
from polyconf.Monomer import Monomer
from polyconf.Polymer import Polymer
from polyconf.PDB import PDB

star = Polymer(Monomer('monomers/extended_neopentane.pdb'))

# join first layer of ethanols
first = star.maxresid()
star.extend(
    Monomer('monomers/extended_ethanol.pdb'),
    n=first,
    nn=star.newresid(),
    names=dict(P='O1', Q='O1', R='C1', S='C1'),
    joins=[('C1', 'O1')],
)
star.extend(
    Monomer('monomers/extended_ethanol.pdb'),
    n=first,
    nn=star.newresid(),
    names=dict(P='O1', Q='O2', R='C1', S='C3'),
    joins=[('C3', 'O1')],
)
star.extend(
    Monomer('monomers/extended_ethanol.pdb'),
    n=first,
    nn=star.newresid(),
    names=dict(P='O1', Q='O3', R='C1', S='C4'),
    joins=[('C4', 'O1')],
)
star.extend(
    Monomer('monomers/extended_ethanol.pdb'),
    n=first,
    nn=star.newresid(),
    names=dict(P='O1', Q='O4', R='C1', S='C5'),
    joins=[('C5', 'O1')],
)
star.renamer(first, 'O1 O2 O3 O4 H12 H1 H6 H9')
for i in range(2, 6):
    star.renamer(i, 'C1')

# add second layer of ethanol
star.extend(
    Monomer('monomers/extended_ethanol.pdb'),
    n=first+1,
    nn=star.newresid(),
    names=dict(P='O1', Q='O2', R='C1', S='C3'),
    joins=[('C3', 'O1')],
)
star.extend(
    Monomer('monomers/extended_ethanol.pdb'),
    n=first+2,
    nn=star.newresid(),
    names=dict(P='O1', Q='O2', R='C1', S='C3'),
    joins=[('C3', 'O1')],
)
star.extend(
    Monomer('monomers/extended_ethanol.pdb'),
    n=first+3,
    nn=star.newresid(),
    names=dict(P='O1', Q='O2', R='C1', S='C3'),
    joins=[('C3', 'O1')],
)
star.extend(
    Monomer('monomers/extended_ethanol.pdb'),
    n=first+4,
    nn=star.newresid(),
    names=dict(P='O1', Q='O2', R='C1', S='C3'),
    joins=[('C3', 'O1')],
)
for i in range(2, 6):
    star.renamer(i, 'O2 H8')
for i in range(6, 10):
    star.renamer(i, 'C1')

# add third layer of ethanol
star.extend(
    Monomer('monomers/extended_ethanol.pdb'),
    n=first+5,
    nn=star.newresid(),
    names=dict(P='O1', Q='O2', R='C1', S='C3'),
    joins=[('C3', 'O1')],
)
star.extend(
    Monomer('monomers/extended_ethanol.pdb'),
    n=first+6,
    nn=star.newresid(),
    names=dict(P='O1', Q='O2', R='C1', S='C3'),
    joins=[('C3', 'O1')],
)
star.extend(
    Monomer('monomers/extended_ethanol.pdb'),
    n=first+7,
    nn=star.newresid(),
    names=dict(P='O1', Q='O2', R='C1', S='C3'),
    joins=[('C3', 'O1')],
)
star.extend(
    Monomer('monomers/extended_ethanol.pdb'),
    n=first+8,
    nn=star.newresid(),
    names=dict(P='O1', Q='O2', R='C1', S='C3'),
    joins=[('C3', 'O1')],
)
for i in range(6, 10):
    star.renamer(i, 'O2 H8')
for i in range(10, 14):
    star.renamer(i, 'C1')

# add caps
star.extend(
    Monomer('monomers/extended_methane.pdb'),
    n=first+9,
    nn=star.newresid(),
    names=dict(P='C1', Q='O2', R='O1', S='C3'),
    joins=[('C3', 'C1')],
)
star.extend(
    Monomer('monomers/extended_methane.pdb'),
    n=first+10,
    nn=star.newresid(),
    names=dict(P='C1', Q='O2', R='O1', S='C3'),
    joins=[('C3', 'C1')],
)
star.extend(
    Monomer('monomers/extended_methane.pdb'),
    n=first+11,
    nn=star.newresid(),
    names=dict(P='C1', Q='O2', R='O1', S='C3'),
    joins=[('C3', 'C1')],
)
star.extend(
    Monomer('monomers/extended_methane.pdb'),
    n=first+12,
    nn=star.newresid(),
    names=dict(P='C1', Q='O2', R='O1', S='C3'),
    joins=[('C3', 'C1')],
)
for i in range(10, 14):
    star.renamer(i, 'O2 H8')
for i in range(14, 18):
    star.renamer(i, 'O1 H4')

Saver = PDB(star)
Saver.cleanup()
Saver.save(dummies='X*',fname=f'polymers/star_polymer') 

# dihedrals = star.gen_pairlist(J='N1',K='C1',first_resid=1,last_resid=127,mult=3,same_res=False) 



# # join first layer of ethanols
# first = star.maxresid()
# star.extend(
#     Monomer('monomers/extended_ethanol.pdb'),
#     n=first,
#     nn=2,
#     names=dict(P='O1', Q='O1', R='C1', S='C1'),
#     joins=[('C1', 'O1')],
# )
# star.extend(
#     Monomer('monomers/extended_ethanol.pdb'),
#     n=first,
#     nn=6,
#     names=dict(P='O1', Q='O2', R='C1', S='C3'),
#     joins=[('C3', 'O1')],
# )
# star.extend(
#     Monomer('monomers/extended_ethanol.pdb'),
#     n=first,
#     nn=10,
#     names=dict(P='O1', Q='O3', R='C1', S='C4'),
#     joins=[('C4', 'O1')],
# )
# star.extend(
#     Monomer('monomers/extended_ethanol.pdb'),
#     n=first,
#     nn=14,
#     names=dict(P='O1', Q='O4', R='C1', S='C5'),
#     joins=[('C5', 'O1')],
# )
# star.renamer(first, 'O1 O2 O3 O4 H12 H1 H6 H9')
# for i in range(2, 6):
#     star.renamer(i, 'C1')

# # add second layer of ethanol
# star.extend(
#     Monomer('monomers/extended_ethanol.pdb'),
#     n=2,
#     nn=3,
#     names=dict(P='O1', Q='O2', R='C1', S='C3'),
#     joins=[('C3', 'O1')],
# )
# star.extend(
#     Monomer('monomers/extended_ethanol.pdb'),
#     n=6,
#     nn=7,
#     names=dict(P='O1', Q='O2', R='C1', S='C3'),
#     joins=[('C3', 'O1')],
# )
# star.extend(
#     Monomer('monomers/extended_ethanol.pdb'),
#     n=10,
#     nn=11,
#     names=dict(P='O1', Q='O2', R='C1', S='C3'),
#     joins=[('C3', 'O1')],
# )
# star.extend(
#     Monomer('monomers/extended_ethanol.pdb'),
#     n=14,
#     nn=15,
#     names=dict(P='O1', Q='O2', R='C1', S='C3'),
#     joins=[('C3', 'O1')],
# )
# for i in range(2, 6):
#     star.renamer(i, 'O2 H8')
# for i in range(6, 10):
#     star.renamer(i, 'C1')

# # add third layer of ethanol
# star.extend(
#     Monomer('monomers/extended_ethanol.pdb'),
#     n=3,
#     nn=4,
#     names=dict(P='O1', Q='O2', R='C1', S='C3'),
#     joins=[('C3', 'O1')],
# )
# star.extend(
#     Monomer('monomers/extended_ethanol.pdb'),
#     n=7,
#     nn=8,
#     names=dict(P='O1', Q='O2', R='C1', S='C3'),
#     joins=[('C3', 'O1')],
# )
# star.extend(
#     Monomer('monomers/extended_ethanol.pdb'),
#     n=11,
#     nn=12,
#     names=dict(P='O1', Q='O2', R='C1', S='C3'),
#     joins=[('C3', 'O1')],
# )
# star.extend(
#     Monomer('monomers/extended_ethanol.pdb'),
#     n=15,
#     nn=16,
#     names=dict(P='O1', Q='O2', R='C1', S='C3'),
#     joins=[('C3', 'O1')],
# )
# for i in range(6, 10):
#     star.renamer(i, 'O2 H8')
# for i in range(10, 14):
#     star.renamer(i, 'C1')

# # add caps
# star.extend(
#     Monomer('monomers/extended_methane.pdb'),
#     n=4,
#     nn=5,
#     names=dict(P='C1', Q='O2', R='O1', S='C3'),
#     joins=[('C3', 'C1')],
# )
# star.extend(
#     Monomer('monomers/extended_methane.pdb'),
#     n=8,
#     nn=9,
#     names=dict(P='C1', Q='O2', R='O1', S='C3'),
#     joins=[('C3', 'C1')],
# )
# star.extend(
#     Monomer('monomers/extended_methane.pdb'),
#     n=12,
#     nn=13,
#     names=dict(P='C1', Q='O2', R='O1', S='C3'),
#     joins=[('C3', 'C1')],
# )
# star.extend(
#     Monomer('monomers/extended_methane.pdb'),
#     n=16,
#     nn=17,
#     names=dict(P='C1', Q='O2', R='O1', S='C3'),
#     joins=[('C3', 'C1')],
# )
# for i in range(10, 14):
#     star.renamer(i, 'O2 H8')
# for i in range(14, 18):
#     star.renamer(i, 'O1 H4')
