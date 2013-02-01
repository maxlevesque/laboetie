subroutine read_elec_slope

  use precision_kinds, only: dp, i2b
  use input, only: input_line
  use system, only: elec_slope_x, elec_slope_y, elec_slope_z

  implicit none
  integer(kind=i2b) :: i, j
  real(kind=dp), dimension(3) :: dummy ! dummy array

  ! give external force a non-physical value
  elec_slope_x = -huge(1.0_dp)
  elec_slope_y = -huge(1.0_dp)
  elec_slope_z = -huge(1.0_dp)

  ! read input file
  j = len ( 'elec_slope' )
  do concurrent( i = 1: size( input_line) )
    if( input_line( i)( 1:j) == 'elec_slope' ) read ( input_line (i) (j+4:j+15) , * ) elec_slope_x, elec_slope_y, elec_slope_z
  end do

  ! has one found these lines
  dummy = (/ elec_slope_x, elec_slope_y, elec_slope_z /)
  if( min(elec_slope_x, elec_slope_y, elec_slope_z)==-huge(1.0_dp)) then
    stop 'in read_elec_slope.f90. elec_slope has not been found correctly in input file.'
  end if

end subroutine read_elec_slope
