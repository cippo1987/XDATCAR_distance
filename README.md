# XDATCAR_distance
This small fortran code calculate atomic distances in not orthogonal cell for a set of structures generated from a MD simulation.

As far as I know there is not smart way to find the nearest atomic image of an atom. All 27 possible images should be considered. 
If the closest vertex is known it is possible to check less cases. 

The code works with the XDATCAR output of VASP and one INPUT file which the following structure
```
Number_structures   Number_of_atomic_species
Number_of_atomic_couples
Atom_A Atom_B #couple 1
Atom_A Atom_B #couple 2
...
Atom_A Atom_B # couple number_of_atomic_couple

```
Atom_A/Atom_B are the number of the atom considered as it appears in the POSCAR.

The code reads a whole step of the simulation, saves it in a variable, then calculates the distance between any of the chosen atoms.

As output are printed the radial distribuction function of the closest atom and of all the equivalent atoms in the structure.
