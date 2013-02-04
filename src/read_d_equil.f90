! read number of iterations one does for equilibration (D_equil)
subroutine read_d_equil

  use precision_kinds, only: i2b
  use input, only: input_line
  use system, only: D_equil

  implicit none
  integer(kind=i2b) :: i, j ! dummy

  D_equil = -huge(1_i2b)
  j = len ( 'D_equil' )
  do concurrent( i = 1: size( input_line) )
    if( input_line( i)( 1:j) == 'D_equil' ) read ( input_line (i) (j+4:j+12) , * ) D_equil
  end do
  if( D_equil == -huge(1_i2b) ) then
    print*, 'D_equil has not been found in input file. check read_d_equil.f90. stop'
    stop
  else if( D_equil <= 0 ) then
    print*, 'D_equil should not be <= 0. check read_d_equil.f90. stop'
    stop
  end if

end subroutine read_d_equil
