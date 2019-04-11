# Compiling
The code should be easy to compile, but beware! -ffixed-line-length-none or -ffree-line-length-n with n large is needed because, well, fortran...


# XDATCAR_distance

This small fortran code calculate atomic distances in not orthogonal cell for a set of structures generated from a MD simulation.
The input has to be named ONE_XDATCAR and it is the simple output of a MD run with VASP (https://www.vasp.at/index.php/about-vasp/59-about-vasp).

This is an input example:
```
Name                                  
           1   
    13.764573    0.000000    0.000000
    -6.862702   11.988418    0.000000
     0.029535    0.001621   14.299894
   Si   O    H    Na   Al
  35  87  30   1   1
Direct configuration=     1
   0.77807603  0.77774692  0.10535127
   0.44009437  0.10733477  0.44003690
   [...]
Direct configuration=     2
```
First line is the name, second line is the scale of the structure. Then you find a 3x3 matrix with the cell vectors.
Line 6 lists the atomic species present. And Line 7 reports the number of atoms per specie.

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
