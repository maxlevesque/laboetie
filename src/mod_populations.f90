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
    ! apply force on all fluid nodes and update populations
    do concurrent (l= lbm%lmin: lbm%lmax)
      where (node%nature == fluid)
        n(:,:,:,l) =  lbm%vel(l)%a0*node%solventDensity &
        + lbm%vel(l)%a1*( lbm%vel(l)%coo(x)*(node%solventFlux(x) + f_ext(x)) &
        + lbm%vel(l)%coo(y)*(node%solventFlux(y) + f_ext(y)) &
        + lbm%vel(l)%coo(z)*(node%solventFlux(z) + f_ext(z)) )
        elsewhere
        n(:,:,:,l) =     lbm%vel(l)%a0*node%solventDensity  &
        + lbm%vel(l)%a1*(  lbm%vel(l)%coo(x)*node%solventFlux(x) &
        + lbm%vel(l)%coo(y)*node%solventFlux(y) &
        + lbm%vel(l)%coo(z)*node%solventFlux(z) )
      end where
    end do

    ! check that no population n < 0
    call check_population (n)
  end subroutine update_populations

  subroutine check_population (arrayin)
    implicit none
    real(dp), dimension(:,:,:,:), intent(in) :: arrayin
    real(dp), parameter :: eps=epsilon(1.0_dp)
    if (any(abs(arrayin) > 1+eps)) stop "Critical. Population n_i(r) > 1 somewhere"
    if (any(arrayin < -eps)) stop 'Critical. Population n_i(r) < 0 somewhere.'
  end subroutine check_population

END MODULE POPULATIONS
