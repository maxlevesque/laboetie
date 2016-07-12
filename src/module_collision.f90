module module_collision

  use precision_kinds, only: dp
  implicit none

  private
  public :: collide, check_population

contains

  subroutine collide(n, density, jx, jy, jz, f_ext_x, f_ext_y, f_ext_z)
    use system, only: fluid, solid, node
    use mod_lbmodel, only: lbm
    use module_input, only: getinput
    implicit none
    real(dp), intent(in) :: density(:,:,:), jx(:,:,:), jy(:,:,:), jz(:,:,:)
    real(dp), intent(in) :: f_ext_x(:,:,:), f_ext_y(:,:,:), f_ext_z(:,:,:)
    real(dp), intent(inout) :: n(:,:,:,:)
    integer :: l, lmin, lmax
    real(dp), allocatable, dimension(:) :: a0, a1, cx, cy, cz
    logical, save :: i_know_the_relaxation_time = .false.
    real(dp), save :: relaxation_time
    real(dp), allocatable :: neq(:,:,:)

    ! With a relaxation time of 1, the population relaxes to Boltzmann distribution
    ! at each timestep.
    ! The kinematic viscosity \nu depends upon the relaxation time (named tau), as :
    ! nu = c_s² * delta_t *(tau/delta_t - 1/2)
    ! Thus if tau=delta_t and c_s²=1/3 for D3Q19 lattice, then nu = 1/6 in lattice units
    if(.not. i_know_the_relaxation_time) then
      relaxation_time = getinput%dp('relaxation_time', defaultvalue=1._dp, assert=">0")
      if( relaxation_time < 0.5_dp) error stop "relaxation_time must be > 0.5"
      i_know_the_relaxation_time = .true.
    end if

    lmin=lbm%lmin
    lmax=lbm%lmax
    allocate( a0(lmin:lmax) )
    allocate( a1(lmin:lmax) )
    allocate( cx(lmin:lmax) )
    allocate( cy(lmin:lmax) )
    allocate( cz(lmin:lmax) )
    a0(lmin:lmax) = lbm%vel(lmin:lmax)%a0
    a1(lmin:lmax) = lbm%vel(lmin:lmax)%a1
    cx(lmin:lmax) = lbm%vel(lmin:lmax)%coo(1)
    cy(lmin:lmax) = lbm%vel(lmin:lmax)%coo(2)
    cz(lmin:lmax) = lbm%vel(lmin:lmax)%coo(3)

    ! First compute the Boltzmann (equilibrium) distribution,
    ! then update the populations according to the relaxation time.
    allocate( neq, mold=jx )

    ! we do the collision even on solid nodes, where density=0
    ! since it does not cost much in cputime and help the compiler to vectorize the (implicit) loops.
    do concurrent( l=lmin:lmax )
      neq(:,:,:) = a0(l)*density + a1(l)*( cx(l)*jx + cy(l)*jy + cz(l)*jz )
      n(:,:,:,l) = (1._dp-1._dp/relaxation_time)*n(:,:,:,l) &
                  + (1._dp/relaxation_time)*neq &
                  + a1(l)*( cx(l)*f_ext_x + cy(l)*f_ext_y + cz(l)*f_ext_z )
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
