subroutine read_lncb_slope

  use precision_kinds, only: dp, i2b
  use input, only: input_line
  use system, only: lncb_slope_x, lncb_slope_y, lncb_slope_z

  implicit none
  integer(kind=i2b) :: i, j
  real(kind=dp), dimension(3) :: dummy ! dummy array

  ! give external force a non-physical value
  lncb_slope_x = -huge(1.0_dp)
  lncb_slope_y = -huge(1.0_dp)
  lncb_slope_z = -huge(1.0_dp)

  ! read input file
  j = len ( 'lncb_slope' )
  do concurrent( i = 1: size( input_line) )
    if( input_line( i)( 1:j) == 'lncb_slope' ) read ( input_line (i) (j+4:j+15) , * ) lncb_slope_x, lncb_slope_y, lncb_slope_z
  end do

  ! has one found these lines
  dummy = (/ lncb_slope_x, lncb_slope_y, lncb_slope_z /)
  if( min(lncb_slope_x, lncb_slope_y, lncb_slope_z)==-huge(1.0_dp)) then
    stop 'in read_lncb_slope.f90. lncb_slope has not been found correctly in input file.'
  end if

end subroutine read_lncb_slope
