! In this subroutine, one prints the first lines of codes to stdout.

subroutine print_header

  implicit none
  character(8)  :: date
  character(10) :: time

  call date_and_time ( DATE=date,TIME=time)

  print*,
  print*,
  print*,date(1:4),'/',date(5:6),'/',date(7:8),', ',time(1:2),'h',time(3:4),'m',time(5:6)
  print*,'===================='
  print*,'Laboetie, the Lattice Boltzmann Electrokinetics '
  print*,'================================================'
  print*,
  print*,

end subroutine print_header
