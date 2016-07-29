! Here tracers are droped in the fluid. They may be charged or not. They evoluate in the
! equilibrated fluid and its solutes. They do not change the potential even if they
! have a charge. The idea is to follow them in order to get the velocity auto correlation
! function, while making them not to change the equilibrium properties of the system.
! Imagine a very small droplet of radioactive particles, so few they do not change
! anything to the system, but numerous enough to be followed and make statistics.

SUBROUTINE drop_tracers(n)

    use precision_kinds, only: dp
    USE system, only: elec_slope
    USE moment_propagation, only: init, propagate, deallocate_propagated_quantity!, print_vacf, integrate_vacf!, not_yet_converged
    use module_input, ONLY: getinput

    IMPLICIT NONE
    real(dp), allocatable, intent(inout) :: n(:,:,:,:) ! xyz,l
    integer :: it, maximum_moment_propagation_steps
    logical :: is_converged

    maximum_moment_propagation_steps = getinput%int("maximum_moment_propagation_steps", 0) ! negative value means make it converge
    IF( maximum_moment_propagation_steps == 0) RETURN

    PRINT*
    PRINT*,'Moment propagation'
    PRINT*,'=================='
    PRINT*,'       step           VACF(x)                   VACF(y)                   VACF(z)'
    PRINT*,'       ----------------------------------------------------------------------------------'

    CALL update_tracer_population ! include elec_slope in population n
    elec_slope = 0.0_dp ! turn the electric field off for the rest of mom_prop (included in n)

    ! add electrostatic potential computed by the SOR routine an external contribution
    ! elec_pot(singlx,ly,lz, ch, phi, phi2, t, t_equil);
    ! call elec_pot
    CALL init ! init moment propagation

    !
    ! Propagation in time
    !
    if( maximum_moment_propagation_steps < 0) maximum_moment_propagation_steps = HUGE(1)
    DO it= 1, maximum_moment_propagation_steps
      !  it=0
      !  do while (not_yet_converged(it))
      !   call elec_pot
      CALL propagate (it,is_converged) ! propagate the propagated quantity
      !    if( modulo(it,50000)==0 ) print_vacf
      !    it = it + 1
      IF( is_converged ) exit
    END DO

    PRINT*,

    !  CALL integrate_vacf ! compute the integral over time of vacf in each direction
    !  CALL print_vacf ! print vacf to file vacf.dat
    CALL deallocate_propagated_quantity

contains


  !
  !
  !
  SUBROUTINE update_tracer_population

    use precision_kinds, only: dp, i2b
    use system, only: f_ext, fluid, elec_slope, supercell, lbm, x, y, z, node
    use module_collision, only: check_population
    use module_input, only: input_dp, input_dp3

    implicit none

    integer(i2b) :: l, ll, lu, lx, ly, lz, i,j,k
    type tracer
      real(dp) :: D ! diffusion coefficient
      real(dp) :: q ! charge
    end type
    type(tracer) :: tr

    ll = lbm%lmin
    lu = lbm%lmax
    lx = supercell%geometry%dimensions%indiceMax(x)
    ly = supercell%geometry%dimensions%indiceMax(y)
    lz = supercell%geometry%dimensions%indiceMax(z)

    deallocate(n) ! the purpose here is to prepare n with the good index arrangement to fit correctly in memory with inner loop over l in moment propagation.
    allocate( n(ll:lu,lx,ly,lz) ,source=0._dp)

    tr%D = input_dp('tracer_Db',0._dp)
    if (tr%D<=epsilon(1._dp)) ERROR STOP 'The diffusion coefficient (tracer_Db in input file) is invalid'

    tr%q = input_dp('tracer_z', 0._dp)
    f_ext(:) = input_dp3('f_ext', [0._dp,0._dp,0._dp] )

    !
    ! Apply force on all fluid nodes and update populations
    !
    do concurrent (l=ll:lu, i=1:lx, j=1:ly, k=1:lz)
      if( node(i,j,k)%nature==fluid ) then
        n(l,i,j,k) = lbm%vel(l)%a0 *node(i,j,k)%solventdensity + lbm%vel(l)%a1 &
          *sum( lbm%vel(l)%coo(:)*(node(i,j,k)%solventFlux(:) +f_ext(:) -node(i,j,k)%solventDensity*tr%q*tr%D*elec_slope(:) ) )
      else
        n(l,i,j,k) = lbm%vel(l)%a0 *node(i,j,k)%solventdensity + lbm%vel(l)%a1 &
          *sum( lbm%vel(l)%coo(:)*node(i,j,k)%solventFlux(:))
      end if
    end do

    !
    ! do concurrent( l= lbm%lmin: lbm%lmax )
    !   where( node%nature ==fluid )
    !     n(l,:,:,:) = lbm%vel(l)%a0*node%solventDensity &
    !                + lbm%vel(l)%a1*(&
    !        lbm%vel(l)%coo(x)*(node%solventFlux(x) + f_ext(x) - node%solventDensity*tr%q *tr%D *elec_slope(x)) &
    !      + lbm%vel(l)%coo(y)*(node%solventFlux(y) + f_ext(y) - node%solventDensity*tr%q *tr%D *elec_slope(y)) &
    !      + lbm%vel(l)%coo(z)*(node%solventFlux(z) + f_ext(z) - node%solventDensity*tr%q *tr%D *elec_slope(z)) )
    !   elsewhere
    !     n(l,:,:,:) = lbm%vel(l)%a0*node%solventDensity + lbm%vel(l)%a1*( &
    !                       lbm%vel(l)%coo(x)*node%solventFlux(x) + &
    !                       lbm%vel(l)%coo(y)*node%solventFlux(y) + &
    !                       lbm%vel(l)%coo(z)*node%solventFlux(z) )
    !   end where
    ! end do

  end SUBROUTINE update_tracer_population

end SUBROUTINE drop_tracers
