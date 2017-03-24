! This is the first part of the LaBo code. Here, no forces are applied on solutes,
! and no flux is allowed. The objective here is to find the initial equilibrium
! distribution of the charges.

SUBROUTINE poisson_nernst_planck

  USE precision_kinds, ONLY: dp
  USE system, ONLY: D_equil, sigma, time, elec_slope, lncb_slope, node, c_plus, c_minus, supercell, phi, & ! ADE: I added phi, c_plus and c_minus
                    fluid
  use module_input, only: getinput
  USE io, ONLY: print_everything_related_to_charge_equil
  use constants, only: x,y,z
  use myallocations


  IMPLICIT NONE

  INTEGER :: timestep, timestepmax, i, j, k, geometrie, lx, ly, lz, capacitor
  LOGICAL :: is_converged = .FALSE.
  real(dp), allocatable, dimension(:,:,:) :: phiTMP ! Ade
  real(dp) :: SF, Alpha, lambda_d

  !print*, 'c_plus = ', abs(c_plus)
  open(314, file='output/c_plus_alongZ.dat')
  open(280, file='output/PHI_PNP.dat')

  !
  ! If sigma is 0, salf free system, and thus no need to go continue here
  !
  sigma = getinput%dp('sigma',0._dp)
  IF ( ABS(sigma) <= EPSILON(1._dp) ) RETURN

  PRINT*
  PRINT*,'Poisson + Nernst-Planck'
  PRINT*,'======================='

  ! read the number of iterations one does for the first step of equilibration (D_iter)

  IF (D_equil<0) STOP "D_equil in input file should be >= 0"

  CALL get_timestepmax (timestepmax)

  ! read electrostatic related stuff
  lncb_slope = getinput%dp3("lncb_slope")
  elec_slope = getinput%dp3("elec_slope")
  

  ! Ade : modification 20/03/2017
  lx = supercell%geometry%dimensions%indiceMax(x)
  ly = supercell%geometry%dimensions%indiceMax(y)
  lz = supercell%geometry%dimensions%indiceMax(z)
  geometrie = getinput%int('geometryLabel',-1)
  IF( .NOT. ALLOCATED(phiTMP) ) CALL allocateReal3D(phiTMP)
  IF( .NOT. ALLOCATED(phi) ) CALL allocateReal3D(phi)
  phi = 0.0_dp ! Ade : initialise phi
  capacitor = getinput%int('capacitor',0) ! 1 = true 0 = false
  if(capacitor.NE.1) then
    if( geometrie == 1 ) then   ! slit pore geometry
        Alpha = getinput%dp('Alpha',0._dp)  ! Attention!!!!!! This is dangerous. We should probably  
                                            ! compute its value in another subroutine
        if(mod(lz,2) == 0) then
           SF = real(lz)/2
        else
           SF = real(lz)/2 + 0.5
        endif
        if( lambda_d == 0.0) then ! no salt added
            do i = 1, lx      
                do j = 1, ly  
                  do k = 1, lz ! first and last node are solid
                    phiTMP(i,j,k) = 2.0*log(cos(Alpha*(k-SF))/(cos(Alpha*(lz-2)/2)) ) ! Potential phi analytical solution for a slit
                    !phiTMP(i,j,k) = 2.0*log(cos(Alpha*(k-SF))) ! Potential phi analytical solution for a slit
                  end do
                end do
            end do
        endif
        ! where(node%nature==solid .and.
        phi = phiTMP 
    endif
  end if

    DO k=1,supercell%geometry%dimensions%indiceMax(z)
      write(280,*) k, SUM(phi(:,:,k)) 
    ENDDO
    close(280)
  ! Ade : end modification 20/03/2017

  ! iterate until charges are equilibrated
    timestep = 0
    DO WHILE ((.NOT. is_converged) .AND. (timestep<=timestepmax))
        CALL backup_phi_c_plus_c_minus ! backup potential and solute concentrations from last step
        CALL sor(timestep) ! TODO    ! compute phi with the Successive Overrelation Routine (SOR)
        CALL just_eq_smolu ! solve smoluchowski (diffusion + electrostatic part) ie not a full smolu
        ! monitor evolution of phi, c_plus, c_minus w.r.t. the backup at the beginning of the iteration
        ! this is done only every 10 loops in order not to waste too much time. arbitrary number.

        IF ((time/=-timestepmax .and. modulo(time,10)==0) .or. sigma==0.0_dp ) THEN
            CALL check_charge_distribution_equilibrium (time, is_converged)
        END IF
    END DO

    DO k=supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
       WRITE(314,*) k, SUM(c_plus(:,:,k))
    ENDDO
    close(314)

    ! Ade : remove comments below. This is a debugging test!!!!!!
    !call charge_test ! check charge conservation ! TODO rename to call check_charge_conservation ! check charge conservation every 1000 steps (arbitrary number)
    if( sigma == 0.0_dp ) return
    call print_everything_related_to_charge_equil

    IF (is_converged) THEN
        PRINT*,'Convergence found at step ',time-1,' after',timestep,' steps'
    ELSE
        STOP 'Equilibrium distribution of salts not found'
    END IF



    CONTAINS

        SUBROUTINE get_timestepmax (a)
            IMPLICIT NONE
            INTEGER :: a
            a = getinput%int("timestepmax_for_PoissonNernstPlanck")
            IF (a<=0) STOP "timestepmax_for_PoissonNernstPlanck must be >=1"
        END SUBROUTINE



end subroutine poisson_nernst_planck
