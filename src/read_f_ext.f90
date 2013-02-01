! Here we read the external force applied to the system (the pressure, for instance)
! as wanted by user in input file

SUBROUTINE READ_F_EXT

  use precision_kinds, only: i2b, dp
  use input, only: input_line
  use system, only: f_ext ! rank 3
  use constants, only: x, y, z

  implicit none
  integer(kind=i2b) :: i, j

  ! give external force a non-physical value
  f_ext = -huge(1.0_dp)

  ! read input file
  j = len ( 'f_ext' )
  do concurrent( i = 1: size( input_line) )
    if( input_line( i)( 1:j) == 'f_ext' ) read ( input_line (i) (j+4:j+50) , * ) f_ext(x), f_ext(y), f_ext(z)
  end do

  ! has one found these lines
  if( minval(f_ext)==-huge(1.0_dp)) then
    stop 'f_ext has not been found correctly in input file. stop in read_f_ext.f90.'
  end if

END SUBROUTINE READ_F_EXT
