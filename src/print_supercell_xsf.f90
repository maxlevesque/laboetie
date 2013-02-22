! Print an XSF file of the supercell for visualisation in VMD for instance.
! Type vmd --xsf output/supercell.xsf to visualise it.
subroutine print_supercell_xsf

use precision_kinds, only: i2b, dp
use system, only: Lx, Ly, Lz, inside, fluid, solid
! Lx size of box in angstroms
! nb_solute_sites nombre de sites (pareil que le nombre de lignes dans format xyz)
! xmol is an array which contains les x de tous les sites, il a donc la taille de x_mol(nb_solute_sites)
! tu peux remplacer integer(kind=i2b) par integer tout court.

implicit none

 integer(kind=i2b) :: i, j, k

 open(5,file='output/supercell.xsf')
 100 format (xA)
 101 format (3(xF10.5))
 102 format (xI5,xI1)
 103 format (xI3,3(xxF10.5))
 write(5,100)'# this is the specification file of the supercell'
 write(5,100)'# lines beginning with # are commented. There cannot be comment lines within the sections'
 write(5,100)'# XSF format specifications can be found on the XCrySDen website http://www.xcrysden.org/doc/XSF.html'
 write(5,100)'# I strongly recommends to read this documentation.'
 write(5,*)
 write(5,100)'# for periodic structures one has to begin with word CRYSTAL'
 write(5,100)'CRYSTAL'
 write(5,100)
 write(5,100)'# Then one needs to specify the lattice vectors'
 write(5,100)'# specification of PRIMVEC (in ANGSTROMS) like:'
 write(5,100)'#         ax, ay, az    (first lattice vector)'
 write(5,100)'#         bx, by, bz    (second lattice vector)'
 write(5,100)'#         cx, cy, cz    (third lattice vector)'
 write(5,100)'# pay attention to vectors as they are written in horizontal way which is quite unusual'
 write(5,100)'# for now only orthorhombic structures are allowed (free norms of lattice vectors, all angles are 90 degrees)'
 write(5,100)'PRIMVEC'
 write(5,101) real(Lx), 0., 0.
 write(5,101) 0., real(Ly), 0.
 write(5,101) 0., 0., real(Lz)
 write(5,*)
 write(5,100)'# Then one needs to specify the atoms belonging to the unit cell. '
 write(5,100)'# First number stands for number of atoms in the primitive cell (2 in this case).'
 write(5,100)'# The second number is always 1 for PRIMCOORD coordinates.'
 write(5,100)'# in angstroms and cartesian coordinates'
 write(5,100)'PRIMCOORD'
 write(5,102) lx*ly*lz, 1
do i = 1, lx
  do j = 1, ly
    do k = 1, lz
        if( inside( i, j, k) == solid ) then
          write(5,103)4, real(i-1,dp), real(j-1,dp), real(k-1,dp)
        else if( inside( i, j, k) == fluid ) then
          write(5,103)1, real(i-1,dp), real(j-1,dp), real(k-1,dp)
        else
          stop 'inside should be solid or fluid only in print_supercell_xsf.f90'
        end if
    end do
  end do
end do

! do i = 1, lx*ly*lz
!        if( inside( i, j, k) == solid ) then
!          write(5,103)4, real(i-1,dp), real(j-1,dp), real(k-1,dp)
!        else if( inside( i, j, k) == fluid ) then
!          write(5,103)1, real(i-1,dp), real(j-1,dp), real(k-1,dp)
!        else
!          stop 'inside should be solid or fluid only in print_supercell_xsf.f90'
!        end if
!!  write(5,103) atomic_nbr(i), x_mol(i), y_mol(i), z_mol(i)
! end do

 close(5)
end subroutine print_supercell_xsf
