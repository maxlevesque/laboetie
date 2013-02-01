! Here we give informations to user

subroutine execution_timer

  use precision_kinds, only: sp

  implicit none
  real(kind=sp), dimension(2) :: tarray
  real(kind=sp) :: time

  call ETIME(tarray,time)

  ! print time informations
  print*,'User time in seconds =',tarray(1)
  print*,'System time in seconds =',tarray(2)
  print*,'Run time since start in seconds =',time

end subroutine execution_timer
