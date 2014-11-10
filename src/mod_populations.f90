MODULE POPULATIONS

  use precision_kinds
  use constants, only: x, y, z
  use system, only: supercell, node
  use mod_lbmodel, only: lbm
  implicit none

  private
  public :: update_populations, check_population

contains

  subroutine update_populations
    use system, only: fluid, f_ext, n, solid
    implicit none
    integer(i2b) :: l
    real(dp) :: a0, a1, c(3)
    ! apply force on all fluid nodes and update populations
    do concurrent (l= lbm%lmin:lbm%lmax)
      a0 = lbm%vel(l)%a0
      a1 = lbm%vel(l)%a1
      c(:) = lbm%vel(l)%coo(:)
      where (node%nature == fluid)
        n(:,:,:,l) =  a0*node%solventDensity &
          + a1*( &
            c(x)*(node%solventFlux(x) + f_ext(x)) &
          + c(y)*(node%solventFlux(y) + f_ext(y)) &
          + c(z)*(node%solventFlux(z) + f_ext(z)) )
      elsewhere
        n(:,:,:,l) =  a0*node%solventDensity &
            + a1*( &
              c(x)*node%solventFlux(x) &
            + c(y)*node%solventFlux(y) &
            + c(z)*node%solventFlux(z)  )
      end where
    end do

    ! check that no population n < 0
    ! call check_population (n)
  end subroutine update_populations

  subroutine check_population (arrayin)
    implicit none
    real(dp), dimension(:,:,:,:), intent(in) :: arrayin
    real(dp), parameter :: eps=epsilon(1.0_dp)
    if (any(abs(arrayin) > 1+eps)) stop "Critical. Population n_i(r) > 1 somewhere"
    if (any(arrayin < -eps)) stop 'Critical. Population n_i(r) < 0 somewhere.'
  end subroutine check_population

END MODULE POPULATIONS
