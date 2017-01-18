! Here tracers are droped in the fluid. They may be charged or not. They evoluate in the
! equilibrated fluid and its solutes. They do not change the potential even if they
! have a charge. The idea is to follow them in order to get the velocity auto correlation
! function, while making them not to change the equilibrium properties of the system.
! Imagine a very small droplet of radioactive particles, so few they do not change
! anything to the system, but numerous enough to be followed and make statistics.

subroutine drop_tracers

  use precision_kinds, only: dp
  USE system, only: tmom, tmax, elec_slope, tic_mp, tac_mp, tic_eq_ads, tac_eq_ads
  USE input, only: input_dp
  USE moment_propagation

  IMPLICIT NONE
  integer :: it
  logical :: is_converged

  CALL update_tracer_population ! include elec_slope in population n
  elec_slope = 0.0_dp ! turn the electric field off for the rest of mom_prop (included in n)

  CALL SYSTEM_CLOCK(tic_eq_ads)
  CALL init_adsorption_equilibration

  DO it= 1, tmax
    CALL Adsorption_equilibration(it, is_converged) !Equilibrate the densities of tracers
    IF( is_converged ) exit
  END DO
  CALL SYSTEM_CLOCK(tac_eq_ads)

  CALL SYSTEM_CLOCK(tic_mp)
    is_converged = .false.
    CALL init_propagation
 
  DO it= 1, tmax
    CALL propagate (it, is_converged) !propagate the quantities
    IF( is_converged ) exit
  END DO
  CALL SYSTEM_CLOCK(tac_mp)
  CALL deallocate_propagated_quantity

  CALL SYSTEM_CLOCK(tac_mp)
contains


  subroutine update_tracer_population
    use precision_kinds, only: dp, i2b
    use system, only: n, f_ext, fluid, elec_slope, lbm, x, y, z, node
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
      where( node%nature ==fluid )
        n(:,:,:,l) = lbm%vel(l)%a0*node%solventDensity &
                   + lbm%vel(l)%a1*(&
           lbm%vel(l)%coo(x)*(node%solventFlux(x) + f_ext(x) - node%solventDensity*tr%q *tr%D *elec_slope(x)) &
         + lbm%vel(l)%coo(y)*(node%solventFlux(y) + f_ext(y) - node%solventDensity*tr%q *tr%D *elec_slope(y)) &
         + lbm%vel(l)%coo(z)*(node%solventFlux(z) + f_ext(z) - node%solventDensity*tr%q *tr%D *elec_slope(z)) )
      elsewhere
        n(:,:,:,l) = lbm%vel(l)%a0*node%solventDensity + lbm%vel(l)%a1*( &
                          lbm%vel(l)%coo(x)*node%solventFlux(x) + &
                          lbm%vel(l)%coo(y)*node%solventFlux(y) + &
                          lbm%vel(l)%coo(z)*node%solventFlux(z) )
      end where
    end do

    call check_population(n)   ! check that no population n < 0
  end subroutine update_tracer_population

end subroutine drop_tracers
