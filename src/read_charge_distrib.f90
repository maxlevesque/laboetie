subroutine read_charge_distrib

  use precision_kinds, only: i2b
  use input, only: input_line
  use system, only: charge_distrib

  implicit none
  character(len=*), parameter :: non_physical_init=':-)'
  integer(kind=i2b) :: i
  integer(kind=i2b), parameter :: j=len('charge_distrib')

  ! init
  charge_distrib = non_physical_init

  ! read
  do concurrent( i = 1: size( input_line) )
    if( input_line( i)( 1:j) == 'charge_distrib' ) read ( input_line (i) (j+4:j+6) , * ) charge_distrib
  end do

  ! test
  if( charge_distrib==non_physical_init ) then
    stop 'in read_charge_distrib.f90. charge_distrib not found in input file.'
  end if

  if( charge_distrib(1:3) /= 'int' .and. charge_distrib(1:3)/='sol') then
    stop 'in read_charge_distrib.f90 charge_distrib can only be int or sol for now'
  end if

end subroutine read_charge_distrib
