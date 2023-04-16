@echo off
setlocal EnableDelayedExpansion

rm concatenated_elements.txt
rm concatenated_topologies.txt
rm concatenated_tests.txt

echo molecule_type.py>> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
type molecule_type.py >> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
echo atoms.py >> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
type atoms.py >> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
echo bonds.py>> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
type bonds.py >> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
echo angles.py>> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
type angles.py >> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
echo dihedrals.py>> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
type dihedrals.py >> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
echo pairs.py>> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
type pairs.py >> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
echo exclusions.py>> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt
type exclusions.py >> concatenated_elements.txt
echo ----------------------------->> concatenated_elements.txt

echo topology.py>> concatenated_topologies.txt
echo ----------------------------->> concatenated_topologies.txt
type topology.py >> concatenated_topologies.txt
echo ----------------------------->> concatenated_topologies.txt
echo monomer.py>> concatenated_topologies.txt
echo ----------------------------->> concatenated_topologies.txt
type monomer.py >> concatenated_topologies.txt
echo ----------------------------->> concatenated_topologies.txt
echo polymer.py>> concatenated_topologies.txt
echo ----------------------------->> concatenated_topologies.txt
type polymer.py >> concatenated_topologies.txt
echo ----------------------------->> concatenated_topologies.txt

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
echo ----------------------------->> concatenated_tests.txt
echo ..\tests\test_dihedrals.py>> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt
type ..\tests\test_dihedrals.py >> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt
echo ..\tests\test_topology.py>> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt
type ..\tests\test_monomer.py >> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt
type ..\tests\test_polymer.py >> concatenated_tests.txt
echo ----------------------------->> concatenated_tests.txt

echo All files have been concatenated