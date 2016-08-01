module module_collision

  use precision_kinds, only: dp
  implicit none

  private
  public :: collide, check_population

contains

  ! for explanation of collision rule, see:
  !     Guo et al. Phys. Rev. E 65, 046308 (2002)
  ! In particular Eqs. 3, 4, 5, 17, 19, 20

  subroutine collide(n, jx, jy, jz, f_ext_x, f_ext_y, f_ext_z)
    use system, only: fluid, solid, node
    use mod_lbmodel, only: lbm
    use module_input, only: getinput
    implicit none
    real(dp), intent(in) :: jx(:,:,:), jy(:,:,:), jz(:,:,:)
    real(dp), intent(in) :: f_ext_x(:,:,:), f_ext_y(:,:,:), f_ext_z(:,:,:)
    real(dp), intent(inout) :: n(:,:,:,:)
    integer :: l, lmin, lmax, nx, ny, nz
    logical, save :: i_know_the_relaxation_time = .false.
    real(dp), save :: relaxation_time
    real(dp), allocatable, dimension(:,:,:), save :: neq, ux, uy, uz, density
    real(dp), parameter :: csq=1.0_dp/3._dp
    logical, save :: i_know_i_want_first_order_only = .false.
    logical, save :: first_order_only = .false.
    integer, save :: collision_order

    ! With a relaxation time of 1, the population relaxes to Boltzmann distribution
    ! at each timestep.
    ! The kinematic viscosity \nu depends upon the relaxation time (named tau), as :
    ! nu = c_s² * delta_t *(tau/delta_t - 1/2)
    ! Thus if tau=delta_t and c_s²=1/3 for D3Q19 lattice, then nu = 1/6 in lattice units
    if(.not. i_know_the_relaxation_time) then
      relaxation_time = getinput%dp('relaxation_time', defaultvalue=1._dp, assert=">0")
      if( relaxation_time < 0.5_dp) error stop "relaxation_time must be > 0.5"
      i_know_the_relaxation_time = .true.
      ! do some initialization like allocations you don't want to do each timestep
      nx = ubound(jx,1)
      ny = ubound(jx,2)
      nz = ubound(jx,3)
      allocate( neq(nx,ny,nz) ,source=0._dp)
      allocate( ux(nx,ny,nz) ,source=0._dp)
      allocate( uy(nx,ny,nz) ,source=0._dp)
      allocate( uz(nx,ny,nz) ,source=0._dp)
      allocate( density(nx,ny,nz), source=0._dp)
    end if

    if(.not. i_know_i_want_first_order_only) then
      first_order_only = getinput%log('first_order_only', defaultvalue=.false.)
      i_know_i_want_first_order_only = .true.
      if( first_order_only ) then
          collision_order = 1
      else
          collision_order = 2
      end if
    end if

    lmin=lbm%lmin
    lmax=lbm%lmax

    ! First compute the Boltzmann (equilibrium) distribution,
    ! then update the populations according to the relaxation time.

    associate( cx=>lbm%cx, cy=>lbm%cy, cz=>lbm%cz, a0=>lbm%a0, a1=>lbm%a1, a2=>lbm%a2 )

    density = sum(n,4)

    where(node%nature==fluid)
      ux=jx/density
      uy=jy/density
      uz=jz/density
    else where
      ux=0
      uy=0
      uz=0
    end where

    select case (collision_order)
    case(2)
      do l=lmin,lmax
        where(node%nature==fluid)

          neq = &
            a0(l)*density &
          + a1(l)*( cx(l)*jx + cy(l)*jy + cz(l)*jz ) &
          + a2(l)*( &
                     jx*ux*(cx(l)**2-csq) + jx*uy*cx(l)*cy(l)    + jx*uz*cx(l)*cz(l) &
                   + jy*ux*cy(l)*cx(l)    + jy*uy*(cy(l)**2-csq) + jy*uz*cy(l)*cz(l) &
                   + jz*ux*cz(l)*cx(l)    + jz*uy*cz(l)*cy(l)    + jz*uz*(cz(l)**2-csq) &
                   )

          n(:,:,:,l) = (1._dp-1._dp/relaxation_time)*n(:,:,:,l) &
          + (1._dp/relaxation_time)*neq &
          + (1._dp-1._dp/(2._dp*relaxation_time))*( &
                a1(l)*( (cx(l)-ux)*f_ext_x + (cy(l)-uy)*f_ext_y + (cz(l)-uz)*f_ext_z ) &
              + 2._dp*a2(l)*( cx(l)*ux + cy(l)*uy + cz(l)*uz )*( cx(l)*f_ext_x + cy(l)*f_ext_y + cz(l)*f_ext_z ))
          ! here we could add the second order also for f_ext
        end where
      end do

    case(1)

      do l=lmin,lmax
        where(node%nature==fluid)

          neq(:,:,:) = a0(l)*density + a1(l)*( cx(l)*jx + cy(l)*jy + cz(l)*jz )
          n(:,:,:,l) = (1._dp-1._dp/relaxation_time)*n(:,:,:,l) &
          + (1._dp/relaxation_time)*neq &
          + (1._dp-1._dp/(2._dp*relaxation_time))*( &
                a1(l)*( (cx(l)-ux)*f_ext_x + (cy(l)-uy)*f_ext_y + (cz(l)-uz)*f_ext_z ) )

        end where
      end do

    end select

    end associate

  end subroutine collide












  subroutine check_population (arrayin)
    implicit none
    real(dp), dimension(:,:,:,:), intent(in) :: arrayin
    real(dp), parameter :: eps=epsilon(1.0_dp)
    if (any(abs(arrayin) > 1+eps)) stop "Critical. Population n_i(r) > 1 somewhere"
    if (any(arrayin < -eps)) stop 'Critical. Population n_i(r) < 0 somewhere.'
  end subroutine check_population

end module module_collision
