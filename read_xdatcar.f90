subroutine calc_vectors(lattice_11, lattice_12, lattice_13, lattice_21, lattice_22, lattice_23, lattice_31, lattice_32, lattice_33, aa, bb, cc, alpha, beta, gamma, calpha, cbeta,cgamma)  
implicit none
real(kind=8), Intent(In)    :: lattice_11, lattice_12, lattice_13,lattice_21, lattice_22, lattice_23, lattice_31, lattice_32, lattice_33
real(kind=8), Intent(Out)   :: aa, bb, cc
real(kind=8), Intent(Out)   :: alpha, beta, gamma, calpha, cbeta, cgamma
real(kind=8)                                                            :: vol
 
 aa = lattice_11
 bb = sqrt(lattice_21*lattice_21  + lattice_22*lattice_22)
 cc = sqrt(lattice_31*lattice_31 + lattice_32*lattice_32+lattice_33*lattice_33)
 vol = (lattice_31*(lattice_12*lattice_23-lattice_13*lattice_22) + lattice_32*(lattice_13*lattice_21-lattice_11*lattice_23) + lattice_33*(lattice_11*lattice_22-lattice_12*lattice_21))
 
alpha = acos((lattice_21*lattice_31 + lattice_22*lattice_32 + lattice_23*lattice_33)/(bb*cc))
beta  = acos((lattice_11*lattice_31 + lattice_12*lattice_32 + lattice_13*lattice_33)/(aa*cc))
gamma = acos((lattice_11*lattice_21 + lattice_12*lattice_22 + lattice_13*lattice_23)/(bb*aa))
 
 calpha = cos(alpha)
 cbeta  = cos(beta)
 cgamma = cos(gamma)
end subroutine calc_vectors
!!!!
!
!!!!
subroutine read_structure(UnitNum,FileName,num_atom, atom_position)  
implicit none
integer, Intent(In)                                                     :: UnitNum
character (len=*), Intent(In)                                           :: FileName
integer(kind=8)                                                         :: i
integer(kind=8), Intent(In)                                             :: num_atom
 
real(kind=8),    DIMENSION(num_atom, 3), Intent(InOut)                  :: atom_position
!write(*,*)UnitNum, FileName
open(unit=UnitNum, file=FileName, status='old', action='read' )
atom_position = -1
read(1,*)! SKIP   !read position of first specie
 
do i=1, num_atom
    read(1,*)atom_position(i,1), atom_position(i,2), atom_position(i,3)    
end do
 
return
end subroutine read_structure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
subroutine rdf(vector_rdf, lung, number_of_couples, couple)
implicit none
integer(kind=8), Intent(InOut)  :: vector_rdf(number_of_couples,600)
real(kind=8), Intent(In)        :: lung
integer(kind=8),Intent(In)      :: number_of_couples, couple
real                    :: bin
integer                 :: position
!write(*,*)lung, "lung"
bin = 0.1
position =  INT(lung/bin) + 1
!write(*,*) "chiamata a rdf: position", position, "bin", bin, "lung", lung
vector_rdf(couple, position)=vector_rdf(couple, position) + 1
RETURN
end subroutine rdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setupvector(vector_rdf,number_of_couples)
implicit none
integer(kind=8), Intent(InOut)      :: vector_rdf(number_of_couples,600)
integer(kind=8),Intent(In)                  :: number_of_couples
vector_rdf=0
RETURN
end subroutine setupvector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program read_xdatcar
implicit none
integer(kind=8)                                         :: num_structu
integer(kind=8)                                         :: num_atom
integer(kind=8)                                         :: num_species, atoma, atomb
integer(kind=8), DIMENSION(:),    ALLOCATABLE           :: atom_a, atom_b
integer(kind=8)                                         :: number_of_couples
 
  real(kind=8),    DIMENSION(:, :), ALLOCATABLE              :: atom_position
integer(kind=8), DIMENSION(:),    ALLOCATABLE              :: element_quantity
integer(kind=8), DIMENSION(3)                             :: per
  real(kind=8)                                          :: lattice_11, lattice_12, lattice_13
real(kind=8)                                            :: lattice_21, lattice_22, lattice_23
real(kind=8)                                            :: lattice_31, lattice_32, lattice_33
  real(kind=8)                                          :: distance, min_distance, bond
  real(kind=8)                                          :: aa,bb,cc
real(kind=8)                                            :: alpha, beta, gamma
real(kind=8)                                            :: calpha, cbeta, cgamma
integer(kind=8)                                         :: i , ii, jj, kk, jkl, write_i
integer(kind=8)                                         :: unit_output
 
integer(kind=8), DIMENSION(:, :), ALLOCATABLE                           :: rdf_minimo, rdf_all, vector_rdf
!integer(kind=8), DIMENSION(600)                            :: rdf_minimo, rdf_all, vector_rdf
 
real(kind=8)                                            :: x_axis
! DEFINE FILES
open (unit =  1, file="ONE_XDATCAR", status='old', action='read')  ! trajectory as in HEADER + ALL THE CONFIGURATIONS
open (unit =  2, file="INPUT", status='old', action='read')  ! trajectory as in HEADER + ALL THE CONFIGURATIONS
 
per(1) = -1
per(2) = 0
per(3) = 1
read(1,*) !num_structu, num_species
read(1,*) !skip
read(1,*) lattice_11, lattice_12, lattice_13
read(1,*) lattice_21, lattice_22, lattice_23
read(1,*) lattice_31, lattice_32, lattice_33
read(1,*) ! skip elements name  
 
!READ INPUT
read(2,*) num_structu, num_species
read(2,*)number_of_couples
 
ALLOCATE(rdf_minimo(number_of_couples,600))
ALLOCATE(rdf_all(number_of_couples,600))
ALLOCATE(vector_rdf(number_of_couples,600))
ALLOCATE(atom_a(number_of_couples))
ALLOCATE(atom_b(number_of_couples))
atom_a = 0
atom_b = 0
 
do jkl=1, number_of_couples
read(2,*)atom_a(jkl),atom_b(jkl)
end do
 
 
 
rdf_minimo = 0
rdf_all = 0
!RESET VECTOR RDF
call setupvector(vector_rdf,number_of_couples)
call calc_vectors(lattice_11, lattice_12, lattice_13, lattice_21, lattice_22, lattice_23, lattice_31, lattice_32, lattice_33, aa, bb, cc, alpha, beta, gamma, calpha, cbeta,cgamma)
! THIS SUB calculates the lenght of the cell and the angles.
!write(*,*)calpha
ALLOCATE(element_quantity(num_species))
element_quantity = 0
SELECT CASE (num_species) ! I don't know how many elements I do have.
   CASE (1)
read(1,*) element_quantity(1)
   CASE (2)
read(1,*) element_quantity(1), element_quantity(2)
   CASE (3)
read(1,*) element_quantity(1), element_quantity(2), element_quantity(3)
   CASE (4)
read(1,*) element_quantity(1), element_quantity(2), element_quantity(3), element_quantity(4)
   CASE (5)
read(1,*) element_quantity(1), element_quantity(2), element_quantity(3), element_quantity(4), element_quantity(5)
   CASE (6)
read(1,*) element_quantity(1), element_quantity(2), element_quantity(3), element_quantity(4), element_quantity(5), element_quantity(6)
   CASE (7)
read(1,*) element_quantity(1), element_quantity(2), element_quantity(3), element_quantity(4), element_quantity(5), element_quantity(6), element_quantity(7)
   CASE (8)
read(1,*) element_quantity(1), element_quantity(2), element_quantity(3), element_quantity(4), element_quantity(5), element_quantity(6), element_quantity(7), element_quantity(8)
END SELECT
!write(*,*)'num species', num_species
num_atom = 0 ! RESET variable
do i=1,num_species
    num_atom = num_atom + element_quantity(i)
end do
!write(*,*)"num atom", num_atom
ALLOCATE(atom_position(num_atom,3) )
atom_position = 0
!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!
!!!    START   !!!
!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!
do i=1,num_structu
 
    call read_structure(1,"ONE_XDATCAR",num_atom, atom_position) ! 1 is the file unit... !ONE_XDATCAR is the file name
    !write(*,*)distance(1,2, aa, bb, cc, calpha, cbeta, cgamma, atom_position, num_atom)
    ! After one structure is read, calculate d, print it.
    do jkl=1,number_of_couples
        atoma = atom_a(jkl)
                atomb = atom_b(jkl)
 
        min_distance = 10000.45   !Necessary to void random number
            do ii=1,3
                do jj=1,3
                    do kk=1,3
                       
                        atom_position(atomb,1) = (atom_position(atomb,1) + per(ii))     ! generate all the 27 possible distance. Not very clever.
                        atom_position(atomb,2) = (atom_position(atomb,2) + per(jj))
                        atom_position(atomb,3) = (atom_position(atomb,3) + per(kk))
                       
                        bond=distance(atoma, atomb, aa, bb, cc, calpha, cbeta, cgamma, atom_position, num_atom)
                        if ( (bond.lt.min_distance) ) then      !.and.(pbclung.gt.0)) then  !write(*,*)"MIN", bond, min_distance
                            min_distance = bond
                        end if  !write(*,*)"this is bond, which is passed as lung", bond
                        !subroutine rdf(vector_rdf, lung, number_of_couples, couple)
                        call rdf(rdf_all, bond, number_of_couples, jkl)                 !write(*,*)"b", bond
                       
                        atom_position(atomb,1) = (atom_position(atomb,1) - per(ii))     ! remove all the 27 possible distance. Not very clever.
                        atom_position(atomb,2) = (atom_position(atomb,2) - per(jj))
                        atom_position(atomb,3) = (atom_position(atomb,3) - per(kk))
 
 
                    end do  !ii
                end do  !! jj
            end do !! kk                !write(*,*)"minima",min_distance
        call rdf(rdf_minimo, min_distance, number_of_couples, jkl)                  !write(*,*)"this is min_distance, which is passed as lung", min_distance
    end do   ! jkl
   
end do ! num_structure
 
    do jkl=1, number_of_couples
        unit_output = 10 + jkl
        x_axis = 0.0
        do write_i=1,600
            write(unit_output,*)x_axis, rdf_minimo(jkl,write_i)
            write(unit_output+30,*)x_axis, rdf_all(jkl,write_i)
            x_axis = x_axis + 0.1
        end do ! jkl writing files
    end do ! i writing files
 
end program read_xdatcar
!!!!!!!!!BEGIN FUNCTION
 
real(kind=8) function distance(atom1, atom2, aa, bb, cc, calpha, cbeta, cgamma, atom_position, num_atom)
implicit none
 
integer(kind=8),      intent(in)                                                :: atom1, atom2, num_atom
real(kind=8), intent(in)                                                        :: aa, bb, cc, calpha, cbeta, cgamma
real(kind=8),    DIMENSION(num_atom,3), Intent(In)                              :: atom_position
real(kind=8)                                                                    :: dfx, dfy, dfz
 dfx = (atom_position(atom1,1) - atom_position(atom2,1))
dfy = (atom_position(atom1,2) - atom_position(atom2,2))
dfz = (atom_position(atom1,3) - atom_position(atom2,3))
 distance = sqrt(aa*aa*dfx*dfx + bb*bb*dfy*dfy + cc*cc*dfz*dfz + 2*bb*cc*calpha*dfy*dfz + 2*cc*aa*cbeta*dfz*dfx + 2*aa*bb*cgamma*dfx*dfy)
 RETURN  !no need of return if contains is used
end function distance
