module module_collision

  use precision_kinds
  implicit none

  private
  public :: collide, check_population

contains

  pure subroutine collide(n, density, jx, jy, jz, f_ext_x, f_ext_y, f_ext_z)
    use system, only: fluid, solid, node
    use mod_lbmodel, only: lbm
    implicit none
    real(dp), intent(in) :: density(:,:,:), jx(:,:,:), jy(:,:,:), jz(:,:,:)
    real(dp), intent(in) :: f_ext_x(:,:,:), f_ext_y(:,:,:), f_ext_z(:,:,:)
    real(dp), intent(inout) :: n(:,:,:,:)
    integer :: l
    real(dp) :: a0, a1, cx, cy, cz

    ! apply force on all fluid nodes and update populations
    do concurrent (l= lbm%lmin:lbm%lmax)
      a0 = lbm%vel(l)%a0
      a1 = lbm%vel(l)%a1
      cx = lbm%vel(l)%coo(1)
      cy = lbm%vel(l)%coo(2)
      cz = lbm%vel(l)%coo(3)
      where (node%nature == fluid)
        n(:,:,:,l) =  a0*density &
          + a1*( &
            cx*(jx + f_ext_x) &
          + cy*(jy + f_ext_y) &
          + cz*(jz + f_ext_z) )
      end where
    end do
    ! call check_population (n) ! check that no population n < 0
  end subroutine collide

  subroutine check_population (arrayin)
    implicit none
    real(dp), dimension(:,:,:,:), intent(in) :: arrayin
    real(dp), parameter :: eps=epsilon(1.0_dp)
    if (any(abs(arrayin) > 1+eps)) stop "Critical. Population n_i(r) > 1 somewhere"
    if (any(arrayin < -eps)) stop 'Critical. Population n_i(r) < 0 somewhere.'
  end subroutine check_population

end module module_collision
