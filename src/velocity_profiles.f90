SUBROUTINE VELOCITY_PROFILES (time)
  use precision_kinds, only : i2b, dp
  use system, only: jx, rho
  use constants, only : x
  implicit none
  integer(i2b), intent(in) :: time

  open(11, file= 'output/velocity_profile')

  write(11,*)
  write(11,*)'# timestep = ',time

  call two_walls_normal_to_z_flux_along_x

  close(11)

  contains


SUBROUTINE two_walls_normal_to_z_flux_along_x
  implicit none
  integer(i2b) :: imin, jmin, k
! velocity profile v_x(z) for
!   z
!   | 
!   |
!   |======================= infinite wall along x
!   |
!   |       flux along x only
!   |
!   |---------------------------> x
!   |
!   |       velocity profile v_y=v_z=0, v_x function of z
!   |
!   |======================= infinite wall along x
!   |
  imin = lbound(jx,1)
  jmin = lbound(jx,2)
  do k = lbound(jx,3), ubound(jx,3)
    write(11,*) k, (jx(imin,jmin,k) )/rho(imin,jmin,k)
  end do

END SUBROUTINE two_walls_normal_to_z_flux_along_x

END SUBROUTINE VELOCITY_PROFILES

