! This subroutine writes an intent(in) array called 'array' to 'filename' in the file format .cube that handles 3-dimensional data.
! The standard .cube file can be read by vmd for instance. It is meant to vizualize periodic supercells.
! see for instance http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html
! or http://paulbourke.net/dataformats/cube/
! It is a very easy format to deal with because it's very natural to write and read.
! Inputs needed : array to write, filename, number of points in each direction, length of each supercell vector, atomic number of the
! atoms and their coordinates. The total number of atoms is determined as the number of different coordinates ie as size(x_mol)
! output : a file called 'filename' containing 3dimensional data in .cube format
! written by Maximilien Levesque, while in post doc at Ecole Normale Superieure, Paris in Daniel Borgis's theoretical chemistry group
! 20110919  Maximilien Levesque, clean version for Virginie M.

subroutine write_to_cube_file ( Lx, Ly, Lz, array, filename )

  use precision_kinds , only : dp, i2b ! defines the precision associated to simple and double

  implicit none
  integer(kind=i2b) :: i, j, k!, n ! dummy variables
  integer(kind=i2b), intent(in) :: lx, ly, lz
  character(*), intent(in) :: filename ! filename of .cube file. For example : "density.cube"
  real(kind=dp), intent(in), dimension ( lx, ly, lz ) :: array ! array printed in .cube file
  real(kind=dp), parameter :: angtobohr = 1.889725989_dp ! 1 Bohr = 1.889725989 Ang. Necessary because of VMD understanding of lengths

  ! define formats for writing to file
!  104 format ( xI3 , xA , 3(xxF10.5) )
  102 format ( 1(xF10.5) )

  ! first open the file you want to print array in. It's a formatted file.
  open(10, file = filename , form = 'formatted' )
  write(10,*) ' CPMD CUBE FILE.' ! default text
  write(10,*) ' OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z' ! default text, for remembering

  ! write the total number of sites and a default text nobody knows it meaning
  write( 10 , * ) 1 ,' 0.0 0.0 0.0 ' ! 0 0 0 ou Lx/2 Ly/2 Lz/2 ?  Size(x_mol) is the total number of atoms

  ! write primary vectors
  write( 10 , * ) lx, lx*angtobohr, ' 0.0 0.0'
  write( 10 , * ) ly, ' 0.0 ' , ly* angtobohr , ' 0.0'
  write( 10 , * ) lz, ' 0.0 0.0 ' , lz* angtobohr

  ! write the atoms and their coordinates in Bohr
  write(10,*) '1 0.0 0.0 0.0 0.0'
!    write( 10 , 104 ) atomic_nbr ( n ) , ' 0.0 ' , x_mol ( n ) * angtobohr , y_mol ( n ) * angtobohr , z_mol ( n ) * angtobohr
!  end do

  ! write .cube file. One value per line. As said earlier, run over x then y then z. The one varying the most rapidly is z.
  do i = 1 , lx
    do j = 1 , ly
      do k = 1 , lz
        write ( 10 , 102 ) array ( i , j , k )
      end do
    end do
  end do

  ! close the .cube file called filename
  close(10)

  ! warn user
  write ( * , * ) filename , ' written'

end subroutine write_to_cube_file
