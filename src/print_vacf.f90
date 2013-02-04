!subroutine print_vacf

!  use system, only: dp, i2b, vacf, tmax, tmom
!  use constants, only: x, y, z

!  implicit none
!  integer(kind=i2b) :: t

!  open(unit=99,file='output/vacf.dat')
!  
!  if( ubound(vacf,2) > tmax-tmom+1) stop 'ubound vacfx too low !?'

!  do t= 0, tmax-tmom
!    write(99,*) t, vacf(x,t), vacf(y,t), vacf(z,t)
!  end do

!  print*,'wrote output/vacf.dat'

!  close(99)

!end subroutine print_vacf
