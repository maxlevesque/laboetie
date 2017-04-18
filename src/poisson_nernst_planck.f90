! This is the first part of the LaBo code. Here, no forces are applied on solutes,
! and no flux is allowed. The objective here is to find the initial equilibrium
! distribution of the charges.

SUBROUTINE poisson_nernst_planck

  USE precision_kinds, ONLY: dp
  USE system, ONLY: D_equil, tot_sol_charge, time, elec_slope, lncb_slope, node, c_plus, c_minus, supercell, phi, & ! ADE: I added phi, c_plus and c_minus
                    fluid, bjl, lambda_D, solid
  use module_input, only: getinput
  USE io, ONLY: print_everything_related_to_charge_equil
  use constants, only: x,y,z, pi
  use myallocations


  IMPLICIT NONE

  INTEGER :: timestep, timestepmax, i, j, k, geometrie, lx, ly, lz, capacitor,&
             slit_sol, m, ENNE
  LOGICAL :: is_converged = .FALSE.
  real(dp), allocatable, dimension(:,:,:) :: phiTMP ! Ade
  real(dp) :: SF, Alpha, FACT1, FACT2, FLAT

  !print*, 'c_plus = ', abs(c_plus)
  open(314, file='output/c_plus_alongZ.dat')
  open(315, file='output/c_minus_alongZ.dat')
  open(316, file='output/c_plus_beforeLoop.dat')
  open(317, file='output/c_plusAlongTime.dat')
  open(280, file='output/PHI_PNP1.dat')
  open(281, file='output/PHI_PNP2.dat')
  open(282, file='output/PHIijk.dat')

  !
  ! If tot_sol_charge is 0, salf free system, and thus no need to go continue here
  !
  tot_sol_charge = getinput%dp('tot_sol_charge',0._dp)
  IF ( ABS(tot_sol_charge) <= EPSILON(1._dp) ) RETURN

  PRINT*
  PRINT*,'Poisson + Nernst-Planck'
  PRINT*,'======================='

  ! read the number of iterations one does for the first step of equilibration (D_iter)
  CALL get_timestepmax (timestepmax)
  D_equil = timestepmax ! Ade: 30/03/2017 D_equil was never given a value. Thus it was always
                        ! equal to zero. This explains why the profiles would never iterate. 

  IF (D_equil<0) STOP "D_equil in input file should be >= 0"
  print*, 'D_equil = ', D_equil


  ! read electrostatic related stuff
  lncb_slope = getinput%dp3("lncb_slope")
  elec_slope = getinput%dp3("elec_slope")
  

  ! Ade : modification 20/03/2017
  lx = supercell%geometry%dimensions%indiceMax(x)
  ly = supercell%geometry%dimensions%indiceMax(y)
  lz = supercell%geometry%dimensions%indiceMax(z)
  geometrie = getinput%int('geometryLabel',-1)
  ! Ade : attention phiTMP here is not the same as in SOR
  IF( .NOT. ALLOCATED(phiTMP) ) CALL allocateReal3D(phiTMP)
  IF( .NOT. ALLOCATED(phi) ) CALL allocateReal3D(phi)
  phi = 0.0_dp ! Ade : initialise phi
  capacitor = getinput%int('capacitor',0) ! 1 = true 0 = false

  slit_sol = getinput%int('slit_sol',0) ! 1= true 0 = false
  ! Ade : m is used for the initial solution input for c_plus        
    if (geometrie==1)  m = 1
    if (geometrie==12) m = 2
    if (geometrie==13) m = 3
    if (geometrie==14) m = 4
  if(capacitor.NE.1 .and. slit_sol==1) then
    if( geometrie == 1 .or. geometrie==12 .or. geometrie==13 .or. geometrie==14) then   ! slit pore geometry
        Alpha = getinput%dp('Alpha',0._dp)  ! Attention!!!!!! This is dangerous. We should probably  
                                            ! compute its value in another subroutine
        ENNE = lz-(2*m) ! Nb of fluid nodes
        if(mod(lz,2) == 0) then
           SF = real(lz)/2
        else
           SF = real(lz)/2 + 0.5
        endif
        if( lambda_d == 0.0) then ! no salt added
            do i = 1, lx      
                do j = 1, ly  
                  do k = m+1, lz-m ! first and last node are solid
                    ! Ade : this is the correct solution. The one underneath was the one in Amael's thesis 
                    ! which unfortunately had a mistake. 
                    phiTMP(i,j,k) = 2.0*log(cos(Alpha*(k-SF))/(cos(Alpha*(real(ENNE)*0.5)) ))! Potential phi analytical solution for a slit
                    !phiTMP(i,j,k) = 2.0*log(cos(Alpha*(k-SF))) ! Potential phi analytical solution for a slit
                  end do
                end do
            end do
            FLAT =  2.0*log(cos(Alpha*(m-SF))/(cos(Alpha*(real(ENNE)*0.5)) ))
            !FLAT =  2.0*log(cos(Alpha*(m-SF)))
        else ! salt added
            FACT1 = 2.0*pi*tot_sol_charge*bjl*lambda_D/(lx*ly)
            FACT2 = ( sinh(real(ENNE)/(2.0*lambda_D)) )
            do i = 1, lx      
                do j = 1, ly  
                  do k = m+1, lz-m ! first and last node are solid
                        phiTMP(i,j,k) = FACT1 * cosh( (k-SF)/lambda_D )/( FACT2 )
                  end do
                end do
            end do
            FLAT =  FACT1 * cosh( (m-SF)/lambda_D )/( FACT2 )
        endif
        phi = phiTMP 
        where(node%nature==solid)
            phi = FLAT
        end where
    endif
  end if

    DO k=1,supercell%geometry%dimensions%indiceMax(z)
      write(280,*) k, SUM(phi(:,:,k)), sum(phiTMP(:,:,k))
    ENDDO
    close(280)
  ! Ade : end modification 20/03/2017
    
    DO k=supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
        WRITE(316,*) k, c_plus(:,:,k)
    ENDDO
    close(316)

  ! iterate until charges are equilibrated
  ! Ade : modifications on 23/03/2017
    time = 0
    DO WHILE ((.NOT. is_converged) .AND. (time<=timestepmax))
        CALL backup_phi_c_plus_c_minus ! backup potential and solute concentrations from last step
        CALL sor ! TODO    ! compute phi with the Successive Overrelation Routine (SOR)
        CALL just_eq_smolu ! solve smoluchowski (diffusion + electrostatic part) ie not a full smolu
        ! monitor evolution of phi, c_plus, c_minus w.r.t. the backup at the beginning of the iteration
        ! this is done only every 10 loops in order not to waste too much time. arbitrary number.
        DO k=supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
        WRITE(317,*) k, SUM(c_plus(:,:,k))
        ENDDO

        IF ((time/= timestepmax .and. modulo(time,10)==0) .or. tot_sol_charge==0.0_dp ) THEN ! Ade : removed minus sign in front of timestepmax
            CALL check_charge_distribution_equilibrium (time, is_converged)
        END IF
        time = time + 1 ! Ade : missing incrementing number
    END DO
    ! ---------------------------------------------------
    ! Ade : POSTPROCESSING AT EQUILIBRIUM 
    ! ---------------------------------------------------
    ! Cations density profile
    DO k=supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
       WRITE(314,*) k, SUM(c_plus(:,:,k))
       WRITE(315,*) k, SUM(c_minus(:,:,k))
    ENDDO
    ! Potential PHI profile
     
    DO i=1,supercell%geometry%dimensions%indiceMax(x)
        DO j=1,supercell%geometry%dimensions%indiceMax(y)
            DO k=1,supercell%geometry%dimensions%indiceMax(z)
                write(282,*) i,j,k, phi(i,j,k) 
            ENDDO
        ENDDO
    ENDDO
    DO k=1,supercell%geometry%dimensions%indiceMax(z)
      write(281,*) k, SUM(phi(:,:,k)), phi(1,1,k)
    ENDDO

    ! Ade : remove comments below. This is a debugging test!!!!!!
    !call charge_test ! check charge conservation ! TODO rename to call check_charge_conservation ! check charge conservation every 1000 steps (arbitrary number)
    if( tot_sol_charge == 0.0_dp ) return
    call print_everything_related_to_charge_equil

    IF (is_converged) THEN
        PRINT*,'Convergence found at step ',time,' after',time,' steps'
    ELSE
        STOP 'Equilibrium distribution of salts not found'
    END IF

    close(314)
    close(315)
    close(317)
    close(281)
    close(282)


    CONTAINS

        SUBROUTINE get_timestepmax (a)
            IMPLICIT NONE
            INTEGER :: a
            a = getinput%int("timestepmax_for_PoissonNernstPlanck",1000000)
            IF (a<=0) STOP "timestepmax_for_PoissonNernstPlanck must be >=1"
        END SUBROUTINE



end subroutine poisson_nernst_planck
