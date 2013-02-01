subroutine read_lambda_d

  use precision_kinds
  use input, only: input_line
  use system, only: lambda_d

  implicit none
  integer(kind=i2b) :: i
  integer(kind=i2b), parameter :: j=len ( 'lambda_D' )
  real(kind=dp), parameter :: non_physical_init= -huge(1.0_dp)

  lambda_D = non_physical_init
  do concurrent( i = 1: size( input_line) )
    if( input_line( i)( 1:j) == 'lambda_D' ) read ( input_line (i) (j+4:j+8) , * ) lambda_d
  end do

  ! test
  if(lambda_d==non_physical_init) then
    stop 'in read_lambda_d.f90. lambda_d not found in input file.'
  end if

end subroutine read_lambda_d
