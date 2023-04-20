@echo off
setlocal EnableDelayedExpansion

rm concatenated_elements.txt
rm concatenated_elements2.txt
rm concatenated_elements3.txt
rm concatenated_topologies.txt
rm concatenated_topologies2.txt
rm concatenated_tests.txt
rm concatenated_tests2.txt
rm concatenated_visualizations.txt


echo molecule_type.py>> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
type molecule_type.py>> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
echo atoms.py >> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
type atoms.py >> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt


echo bonds.py>> concatenated_elements3.txt
echo ----------------------------->> concatenated_elements3.txt
type bonds.py >> concatenated_elements3.txt
echo ----------------------------->> concatenated_elements3.txt
echo angles.py>> concatenated_elements3.txt
echo ----------------------------->> concatenated_elements3.txt
type angles.py >> concatenated_elements3.txt
echo ----------------------------->> concatenated_elements3.txt


echo dihedrals.py>> concatenated_elements2.txt
echo ----------------------------->> concatenated_elements2.txt
type dihedrals.py >> concatenated_elements2.txt
echo ----------------------------->> concatenated_elements2.txt
echo pairs.py>> concatenated_elements2.txt
echo ----------------------------->> concatenated_elements2.txt
type pairs.py >> concatenated_elements2.txt
echo ----------------------------->> concatenated_elements2.txt
echo exclusions.py>> concatenated_elements2.txt
echo ----------------------------->> concatenated_elements2.txt
type exclusions.py >> concatenated_elements2.txt
echo ----------------------------->> concatenated_elements2.txt

echo topology.py>> concatenated_topologies.txt
echo ----------------------------->> concatenated_topologies.txt
type topology.py >> concatenated_topologies.txt

echo ----------------------------->> concatenated_topologies2.txt
echo monomer.py>> concatenated_topologies2.txt
echo ----------------------------->> concatenated_topologies2.txt
type monomer.py >> concatenated_topologies2.txt
echo ----------------------------->> concatenated_topologies2.txt
echo polymer.py>> concatenated_topologies2.txt
echo ----------------------------->> concatenated_topologies2.txt
type polymer.py >> concatenated_topologies2.txt
echo ----------------------------->> concatenated_topologies2.txt

type visualize.py >> concatenated_visualizations.txt
echo ----------------------------->> concatenated_visualizations.txt

echo ..\tests\test_atoms.py>> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt
type ..\tests\test_atoms.py >> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt
echo ..\tests\test_bonds.py>> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt
type ..\tests\test_bonds.py >> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt
echo ..\tests\test_angles.py>> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt
type ..\tests\test_angles.py >> concatenated_tests.txt
echo ----------------------------->> concatenated_tests2.txt
echo ..\tests\test_dihedrals.py>> concatenated_tests2.txt
echo ----------------------------->> concatenated_tests2.txt
type ..\tests\test_dihedrals.py >> concatenated_tests2.txt
echo ----------------------------->> concatenated_tests.txt
echo ..\tests\test_topology.py>> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt
type ..\tests\test_monomer.py >> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt
type ..\tests\test_polymer.py >> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt

echo All files have been concatenated