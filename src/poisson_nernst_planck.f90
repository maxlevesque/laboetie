! This is the first part of the LaBo code. Here, no forces are applied on solutes,
! and no flux is allowed. The objective here is to find the initial equilibrium
! distribution of the charges.

subroutine poisson_nernst_planck

  use precision_kinds, only: i2b, dp
  use system, only: D_equil, sigma, time
  use input, only: input_dp
  implicit none
  logical :: is_converged = .false.

  ! test if salt free system, no solute to equilibrate: return to main
  sigma = input_dp('sigma')
  if( sigma == 0.0_dp ) then ! sigma has already been read in charge initialisation
    print*,'salt free system.'
    ! it remains important to go through all next steps in order to insure the reading of everything, etc.
  end if

  ! read the number of iterations one does for the first step of equilibration (D_iter)
  call read_D_equil

  ! read electrostatic related stuff
  call read_lncb_slope
  call read_elec_slope

  ! iterate until charges are equilibrated
  equiloop: do time = -D_equil, 0 ! time 0 being first step with flux

    ! print time
!    print*,;print*,'time =',time

    ! backup potential and solute concentrations from last step
    call backup_phi_c_plus_c_minus

    ! compute phi with the Successive Overrelation Routine (SOR)
    call sor ! TODO

    ! solve smoluchowski (diffusion + electrostatic part) ie not a full smolu
    call just_eq_smolu ! called just_equ_smolu in c code

    ! check charge conservation every 1000 steps (arbitrary number)
    if(modulo(time,1000)==0) call charge_test ! TODO rename to call check_charge_conservation

    ! monitor evolution of phi, c_plus, c_minus w.r.t. the backup at the beginning of the iteration
    ! this is done only every 10 loops in order not to waste too much time. arbitrary number.
    if((time/=-D_equil .and. modulo(time,10)==0) .or. sigma==0.0_dp ) then
      call check_charge_distribution_equilibrium (time, is_converged)
    end if

    ! if equilibrium is found, exit
    if(is_converged) exit equiloop ! for now always .false. but that's a simple and efficient criteria to implement in a do while equiloop

  end do equiloop

  ! check charge conservation
  call charge_test

  if( sigma == 0.0_dp ) return

  call print_everything_related_to_charge_equil

  if(is_converged) then
    print*,'convergence found at step ',time-1,' after',abs(-D_equil-time),' steps'
  else
    stop 'at the end of poisson_nernst_planck.f90, charges are still not converged'
  end if


end subroutine poisson_nernst_planck
