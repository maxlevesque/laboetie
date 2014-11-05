! Here tracers are droped in the fluid. They may be charged or not. They evoluate in the
! equilibrated fluid and its solutes. They do not change the potential even if they
! have a charge. The idea is to follow them in order to get the velocity auto correlation
! function, while making them not to change the equilibrium properties of the system.
! Imagine a very small droplet of radioactive particles, so few they do not change
! anything to the system, but numerous enough to be followed and make statistics.

subroutine drop_tracers

  use precision_kinds, only: dp
  USE system, only: tmom, tmax, elec_slope
  USE moment_propagation, only: init, propagate, deallocate_propagated_quantity!, print_vacf, integrate_vacf!, not_yet_converged
  USE input, only: input_dp

  IMPLICIT NONE
  integer :: it
  logical :: is_converged

  PRINT*,'       step           VACF(x)                   VACF(y)                   VACF(z)'
  PRINT*,'       ----------------------------------------------------------------------------------'

  CALL update_tracer_population ! include elec_slope in population n
  elec_slope = 0.0_dp ! turn the electric field off for the rest of mom_prop (included in n)

  ! add electrostatic potential computed by the SOR routine an external contribution
  ! elec_pot(singlx,ly,lz, ch, phi, phi2, t, t_equil);
  ! call elec_pot
  CALL init ! init moment propagation

  ! propagate in time
  DO it= 1, tmax-tmom
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


  subroutine update_tracer_population
    use precision_kinds, only: dp, i2b
    use system, only: n, f_ext, fluid, elec_slope, supercell, lbm, x, y, z
    use populations, only: check_population
    use input, only: input_dp
    implicit none
    integer(i2b) :: l
    type tracer
      real(dp) :: D ! diffusion coefficient
      real(dp) :: q ! charge
    end type
    type(tracer) :: tr

    tr%D = input_dp('tracer_Db')
    if (tr%D<=epsilon(1._dp)) stop 'D_tracer, ie tracer_Db in input is invalid'

    tr%q = input_dp('tracer_z')

    ! apply force on all fluid nodes and update populations
    do concurrent( l= lbm%lmin: lbm%lmax )
      where( supercell%node%nature ==fluid )
        n(:,:,:,l) = lbm%vel(l)%a0*supercell%node%solventDensity &
                   + lbm%vel(l)%a1*(&
           lbm%vel(l)%coo(x)*(supercell%node%solventFlux(x) + f_ext(x) &
                              - supercell%node%solventDensity*tr%q *tr%D *elec_slope(x)) &
         + lbm%vel(l)%coo(y)*(supercell%node%solventFlux(y) + f_ext(y) &
                              - supercell%node%solventDensity*tr%q *tr%D *elec_slope(y)) &
         + lbm%vel(l)%coo(z)*(supercell%node%solventFlux(z) + f_ext(z) &
                              - supercell%node%solventDensity*tr%q *tr%D *elec_slope(z)) )
      elsewhere
        n(:,:,:,l) = lbm%vel(l)%a0*supercell%node%solventDensity + lbm%vel(l)%a1*( &
                          lbm%vel(l)%coo(x)*supercell%node%solventFlux(x) + &
                          lbm%vel(l)%coo(y)*supercell%node%solventFlux(y) + &
                          lbm%vel(l)%coo(z)*supercell%node%solventFlux(z) )
      end where
    end do

    call check_population(n)   ! check that no population n < 0
  end subroutine update_tracer_population

end subroutine drop_tracers
