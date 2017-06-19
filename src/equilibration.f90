SUBROUTINE equilibration

    USE precision_kinds, only: i2b, dp, sp
    USE system, only: fluid, supercell, node, lbm, n, pbc, solute_force, t_equil, c_plus, c_minus, phi, Phi_tot
    use module_collision, only: collide
    use module_input, only: getinput
    USE constants, only: x, y, z, zerodp
    USE mod_time, only: tick, tock
    USE myallocations

    implicit none
    integer :: t,i,j,k,l,ip,jp,kp, lmin, lmax, timer(100), g, ng, pdr, pd, ios, px, py, pz, pCoord(3)
    integer :: fluid_nodes, print_frequency, supercellgeometrylabel, tfext, print_files_frequency, GL, print_every
    integer(kind(fluid)), allocatable, dimension(:,:,:) :: nature
    real(dp) :: n_loc, f_ext_loc(3), l2err, target_error, djx, djy, djz, Jxx, Jyy, Jzz
    REAL(dp) :: vmaxx, vmaxy, vmaxz, vmax
    REAL(dp) :: vmaxx_old, vmaxy_old, vmaxz_old, vmax_old
    real(dp), allocatable, dimension(:,:,:) :: density, jx, jy, jz, n_old, jx_old, jy_old, jz_old, f_ext_x, f_ext_y, f_ext_z,& 
                                               F1, F2, F3
    real(dp), allocatable, dimension(:) :: a0, a1
    integer, allocatable, dimension(:) :: cx, cy, cz
    logical :: convergence_reached, compensate_f_ext, convergence_reached_without_fext, convergence_reached_with_fext, err
    REAL(dp), PARAMETER :: eps=EPSILON(1._dp)
    LOGICAL :: write_total_mass_flux
    integer, allocatable :: il(:,:), jl(:,:), kl(:,:), l_inv(:)
    integer(i2b) :: lx, ly, lz, half
    integer(i2b) :: n1,n2,n3 ! Ade : dummy for reading purposes
    character*200 :: ifile, ifile2
    LOGICAL :: RestartPNP = .TRUE.

    open(316, file = "output/soluteForceEqX.dat")
    open(323, file = "output/soluteForceEqY.dat")
    open(324, file = "output/soluteForceEqZ.dat")
    open(1316, file = "output/SFX.dat")
    open(1323, file = "output/SFY.dat")
    open(1324, file = "output/SFZ.dat")
    open(1325, file = "output/SFXtime.dat")
    open(1326, file = "output/SFYtime.dat")
    open(1327, file = "output/SFZtime.dat")
    ! For debugging
    open(1328, file = "output/Ligne2Courant.dat")
    ! end debugging
    open(325, file = "output/phi.dat")
    open(387, file = "output/phiAVANT.dat")
    open(388, file = "output/phiAPRES.dat")
    open(389, file = "output/c_plusTimeEquil.dat")
    open(390, file = "output/PHITimeEquil.dat")
    !------------------------------------------------
    ! Ade : For restart purposes
    !open(1389, file = "output/PHIijk.dat")
    !open(1390, file = "output/C_PLUSijk.dat")
    !open(1391, file = "output/C_MINUSijk.dat")
    !------------------------------------------------

    lmin = lbm%lmin
    lmax = lbm%lmax

    supercellgeometrylabel = supercell%geometry%label ! -1 for solid free cell

    lx = supercell%geometry%dimensions%indiceMax(x)
    ly = supercell%geometry%dimensions%indiceMax(y)
    lz = supercell%geometry%dimensions%indiceMax(z)
    !******************************
    !******* Restart-PNP **********
    !******************************
    RestartPNP = getinput%log("RestartPNP", .TRUE.)
    !if(.not.RestartPNP) then
    !    IF( .NOT. ALLOCATED(phi) ) CALL allocateReal3D(phi)
    !   IF( .NOT. ALLOCATED(c_plus) ) CALL allocateReal3D(c_plus)
    !    IF( .NOT. ALLOCATED(c_minus) ) CALL allocateReal3D(c_minus)
    !    DO i=1,lx
    !     DO j=1,ly
    !      DO k=1,lz
    !        read(1389,*) phi(i,j,k)
    !        read(1390,*) c_plus(i,j,k)
    !        read(1391,*) c_minus(i,j,k)
    !      ENDDO
    !     ENDDO
    !    ENDDO
    !endif
    !******************************
    !******* End of Restart-PNP ***
    !******************************

    !
    ! Print info to terminal every that number of steps
    !
    print_frequency = getinput%int('print_frequency', defaultvalue=max(int(50000/(lx*ly*lz)),1), assert=">0" ) ! this number is my own optimal. To be generalized on strong criteria some day.

    !
    ! WRITE velocity profiles to terminal every that number of steps
    !
    print_files_frequency = getinput%int("print_files_frequency", HUGE(1))

    print_every = getinput%int("print_every", defaultvalue=0) ! reads from lb.in file
                                                                 ! the frequency of printing time
                                                                 ! the default value needs to be changed eventually
    fluid_nodes = count( node%nature==fluid )

    ! Max : I had 1.D-8 before ADE's modification (June 21)
    target_error = getinput%dp("target_error", 1.D-10)

    allocate( density(lx,ly,lz), source=node%solventdensity, stat=ios)
    if (ios /= 0) stop "density: Allocation request denied"

    allocate( l_inv(lmin:lmax) , stat=ios)
    if (ios /= 0) stop "l_inv: Allocation request denied"
    do l = lmin, lmax
        l_inv(l) = lbm%vel(l)%inv
    end do

    !--------------------------------- ADE -----------------------------------------------------------------
    ! ADE : we allocate solute_force
    if(.not.allocated(solute_force)) allocate(solute_force(lx,ly,lz,x:z),source=0.0_dp)
    !--------------------------------- ADE -----------------------------------------------------------------

    allocate( jx     (lx,ly,lz), source=node%solventflux(x))
    allocate( jx_old (lx,ly,lz) )
    allocate( jy     (lx,ly,lz), source=node%solventflux(y))
    allocate( jy_old (lx,ly,lz) )
    allocate( jz     (lx,ly,lz), source=node%solventflux(z))
    allocate( jz_old (lx,ly,lz) )
    jx = 0
    jy = 0
    jz = 0
    jx_old = 0
    jy_old = 0
    jz_old = 0
    djx = 0.0
    djy = 0.0
    djz = 0.0

    OPEN(66, FILE="output/mass-flux_profile_along_z.dat")
    OPEN(67, FILE="output/mass-flux_profile_along_y.dat")
    OPEN(68, FILE="output/mass-flux_profile_along_x.dat")
    WRITE(66,*) "# z, <ρ.v_x>_{x,y}, <ρ.v_y>_{x,y}, <ρ.v_z>_{x,y}"
    WRITE(67,*) "# y, <ρ.v_x>_{x,z}, <ρ.v_y>_{x,z}, <ρ.v_z>_{x,z}"
    WRITE(68,*) "# x, <ρ.v_x>_{y,z}, <ρ.v_y>_{y,z}, <ρ.v_z>_{y,z}"

    OPEN(56, FILE="output/mean-density_profile_along_z.dat")
    OPEN(57, FILE="output/mean-density_profile_along_y.dat")
    OPEN(58, FILE="output/mean-density_profile_along_x.dat")

    allocate( nature (lx,ly,lz), source=node%nature)
    allocate( f_ext_x(lx,ly,lz), source=zerodp)
    allocate( f_ext_y(lx,ly,lz), source=zerodp)
    allocate( f_ext_z(lx,ly,lz), source=zerodp)
    !--------------- Ade ----------------------
    allocate( F1(lx,ly,lz), source=zerodp)
    allocate( F2(lx,ly,lz), source=zerodp)
    allocate( F3(lx,ly,lz), source=zerodp)
    !--------------- Ade ----------------------

    f_ext_loc = zerodp ! this is important and spagetty like... please read carefully before modifying this line
    allocate( cx(lmax), source=lbm%vel(:)%coo(1))
    allocate( cy(lmax), source=lbm%vel(:)%coo(2))
    allocate( cz(lmax), source=lbm%vel(:)%coo(3))
    allocate( a0(lmax), source=lbm%vel(:)%a0)
    allocate( a1(lmax), source=lbm%vel(:)%a1)

    !
    ! Tabulate the index of the node one finishes if one starts from a node and a velocity index l
    ! per direction
    !
    allocate( il(lbm%lmin:lbm%lmax, 1:lx), stat=ios)
    if (ios /= 0) stop "il: Allocation request denied"
    allocate( jl(lbm%lmin:lbm%lmax, 1:ly), stat=ios)
    if (ios /= 0) stop "jl: Allocation request denied"
    allocate( kl(lbm%lmin:lbm%lmax, 1:lz), stat=ios)
    if (ios /= 0) stop "kl: Allocation request denied"
    do l= lmin, lmax
        il(l,:) = [( pbc(i+cx(l),x) ,i=1,lx )]
        jl(l,:) = [( pbc(j+cy(l),y) ,j=1,ly )]
        kl(l,:) = [( pbc(k+cz(l),z) ,k=1,lz )]
    end do

    convergence_reached_without_fext = .false.
    convergence_reached_with_fext = .false.
    compensate_f_ext = getinput%log("compensate_f_ext",.false.)
    if(compensate_f_ext) then
        open(79,file="./output/v_centralnode.dat")
        open(80,file="./output/rho_centralnode.dat")
    endif

    write_total_mass_flux = getinput%log("write_total_mass_flux", .FALSE.)
    IF( write_total_mass_flux ) THEN
        OPEN( 65, FILE="output/total_mass_flux.dat" )
    END IF



    PRINT*
    PRINT*,'Lattice Boltzmann'
    PRINT*,'================='
    PRINT*,'       step'
    PRINT*,'       ----'


    ! ADE : We initialise tfext, which is the time from when the force f_ext is applied
    ! upon a certain number of nodes. 
    tfext = HUGE(tfext)
    ! We also initialise l2err, the convergence error
    l2err = -1

    !
    ! TIME STEPS (in lb units)
    !
    do t=1,HUGE(t)


        !--------------------------------- Ade --------------------------------------------------------------
        if (t<t_equil) then
            f_ext_x = zerodp
            f_ext_y = zerodp
            f_ext_z = zerodp
        else
            f_ext_x = f_ext_loc(1)
            f_ext_y = f_ext_loc(2)
            f_ext_z = f_ext_loc(3)
        endif
        !--------------------------------- Ade --------------------------------------------------------------

        !-------------------------------- Ade -----------------------------------------------------------
        !print*, 2, SUM(c_plus(:,:,2))
        ! This is a debugging test
        !-------------------------------- Ade -----------------------------------------------------------

        !
        ! Print sdtout timestep, etc
        !
        IF( MODULO(t, print_frequency) == 0) PRINT*, t, real(l2err),"(target",real(target_error,4),")"

        !
        ! WRITE velocity profiles
        !
        IF( MODULO(t, print_files_frequency) == 0 .OR. t==1) THEN
            WRITE(66,*)"# timestep",t
            WRITE(67,*)"# timestep",t
            WRITE(68,*)"# timestep",t
            WRITE(56,*)"# timestep",t
            WRITE(57,*)"# timestep",t
            WRITE(58,*)"# timestep",t
            DO k=1,lz
                WRITE(56,*) k, SUM(density(:,:,k))/ MAX( COUNT(density(:,:,k)>eps)  ,1)
                WRITE(66,*) k, SUM(jx(:,:,k)), SUM(jy(:,:,k)), SUM(jz(:,:,k))
            END DO
            DO k=1,ly
                WRITE(57,*) k, SUM(density(:,k,:))/ MAX( COUNT(density(:,k,:)>eps)  ,1)
                WRITE(67,*) k, SUM(jx(:,k,:)), SUM(jy(:,k,:)), SUM(jz(:,k,:))
            END DO
            DO k=1,lx
                WRITE(58,*) k, SUM(density(k,:,:))/ MAX( COUNT(density(k,:,:)>eps)  ,1)
                WRITE(68,*) k, SUM(jx(k,:,:)), SUM(jy(k,:,:)), SUM(jz(k,:,:))
            END DO
            WRITE(66,*)
            WRITE(67,*)
            WRITE(68,*)
        END IF

        !
        ! Compensate_f_ext is an option to have a local force one could apply on a given node (only the central node for now)
        ! compensated by a continuum background force, a little bit like a compensating electric field in a charged supercell.
        !
        ! Ade: By doing so, we can apply a force fx or fy in order to analyse and observe the velocity streamlines around
        ! the particle in the output file v_centralnode.dat.
        !
        if( compensate_f_ext .and. convergence_reached_without_fext) then
            if(px<=0 .or. py<=0 .or. pz<=0) error stop "px, py or pz is not valid in equilibration.f90"
            write(79,*)t-tfext, jx(px,py,pz), jy(px,py,pz), jz(px,py,pz)
            write(80,*)t-tfext, density(px,py,pz)
        end if


        !--------------------------------- ADE -----------------------------------------------------------------
        ! Ade : BEGIN OF MODIF 17/01/17

        ! f_ext is obtained reading input file lb.in (=> pressure gradient)
        ! solute_force is computed in smolu.f90

        write(316,*) 't =', t
        write(323,*) 't =', t
        write(324,*) 't =', t
        DO k=1,lz
            write(316,*) k, SUM(solute_force(:,:,k,1)) ! Ade : The fluid is moving in the y-direction whenever a slit 
                                                        ! case is imposed, as the walls are located at z = 0 and z = L
                                                        ! which is the reason why we are observing F_y(z). 2=>y and k=>z
            write(323,*) k, SUM(solute_force(:,:,k,2)) ! Ade : The fluid is moving in the y-direction whenever a slit 
            write(324,*) k, SUM(solute_force(:,:,k,3)) ! Ade : The fluid is moving in the y-direction whenever a slit 
        END DO


        F1(:,:,:)  = f_ext_x(:,:,:) + solute_force(:,:,:,1)
        F2(:,:,:)  = f_ext_y(:,:,:) + solute_force(:,:,:,2)
        F3(:,:,:)  = f_ext_z(:,:,:) + solute_force(:,:,:,3) 


        !##################
        !# Collision step #
        !##################

        call collide(n, density, jx, jy, jz, F1, F2, F3)
        !call collide(n, density, jx, jy, jz, f_ext_x, f_ext_y, f_ext_z) ! Ade : There is a problem in solute_force

        !--------------------------------- ADE -----------------------------------------------------------------
        
        ! print velocity profile if you need/want it
        ! if( modulo(t, print_frequency) == 0) then
        !    call velocity_profiles(t) ! print velocity profiles
        ! end if

        !
        ! Bounce back (boundpm) to simplify propagation step
        !
        do concurrent(l=lmin:lmax:2)
            do concurrent(k=1:lz)
                kp = kl(l,k)
                !kp=pbc(k+cz(l),z)
                do concurrent(j=1:ly)
                    jp = jl(l,j)
                    !jp=pbc(j+cy(l),y)
                    do concurrent(i=1:lx)
                        ip = il(l,i)
                        !ip=pbc(i+cx(l),x)
                        if( nature(i,j,k) /= nature(ip,jp,kp) ) then
                            n_loc = n(i,j,k,l)
                            n(i,j,k,l) = n(ip,jp,kp,l_inv(l))
                            n(ip,jp,kp,l_inv(l)) = n_loc
                        end if
                    end do
                end do
            end do
        end do

        !###############
        !# PROPAGATION #
        !###############
        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(n,lz,ly,lx,lmin,lmax,cz,cy,cx,il,jl,kl) &
        !$OMP PRIVATE(l,k,j,i,ip,jp,kp,n_old)
        do l=lmin,lmax
            n_old = n(:,:,:,l)
            do k=1,lz
                kp = kl(l,k)
                do j=1,ly
                    jp = jl(l,j)
                    do i=1,lx
                        ip = il(l,i)
                        n(ip,jp,kp,l) = n_old(i,j,k)
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        !
        ! The populations are probabilities and thus must never be negative
        !
        IF( ANY(n<0) ) ERROR STOP "In equilibration, the population n(x,y,z,vel) < 0"

        !
        ! Update densities after the propagation and check it
        ! Densities are the sum of all velocities of a local population
        !
        density = SUM(n,4)

        !###########################
        !# WRITE the total density #
        !###########################
        IF( write_total_mass_flux ) THEN
            WRITE(65,*) t, REAL([  SUM(jx), SUM(jy), SUM(jz)  ])
        END IF

        !
        ! backup moment density (velocities) to test convergence at the end of the timestep
        !
        jx_old = jx
        jy_old = jy
        jz_old = jz

        !IF( MODULO(t, print_frequency) == 0) PRINT*,  t, "before", jz(1,1,1)

        ! update momentum densities after the propagation
        ! this is completely local in space and my be parallelized very well
        ! !$OMP PARALLEL DO DEFAULT(NONE)&
        ! !$OMP PRIVATE(l)&
        ! !$OMP SHARED(lmin,lmax,n,cx,cy,cz)&
        ! !$OMP REDUCTION(+:jx)&
        ! !$OMP REDUCTION(+:jy)&
        ! !$OMP REDUCTION(+:jz)
        ! do l=lmin,lmax
        !     jx = jx +n(:,:,:,l)*cx(l)
        !     jy = jy +n(:,:,:,l)*cy(l)
        !     jz = jz +n(:,:,:,l)*cz(l)
        ! end do
        ! !$OMP END PARALLEL DO
        ! jx=jx/2
        ! jy=jy/2
        ! jz=jz/2
        ! BEN+MAX: 12/07/2016 change the way we integrate n_l*c_l
        !jx=0
        !jy=0
        !jz=0
        !jx=f_ext_x/2._dp
        !jy=f_ext_y/2._dp
        !jz=f_ext_z/2._dp

        jx = F1/2._dp ! Ade : 30/05/17 there was a mistake here. Only f_ext was divided
        jy = F2/2._dp ! by 2.0_dp
        jz = F3/2._dp
        do l=lmin,lmax
            jx = jx +n(:,:,:,l)*cx(l)
            jy = jy +n(:,:,:,l)*cy(l)
            jz = jz +n(:,:,:,l)*cz(l)
        end do

        !IF( MODULO(t, print_frequency) == 0) PRINT*,  t, "after ", jz(1,1,1)
        
        !--------------------------------- ADE -----------------------------------------------------------------

        ! Ade: we need to assign the correct table to node%solventflux as it is being called by other following
        ! subroutines (e.g. advect )
        node%solventdensity = density
        node%solventflux(x) = jx
        node%solventflux(y) = jy
        node%solventflux(z) = jz
        !node%nature = nature ! Ade : Do I actually need this line here ????

        call advect
        call sor ! TODO    ! compute phi with the Successive Overrelation Routine (SOR)
        ! --------------------------- Ade : 13/02/2017 ---------------------------------------------------------------
        !write(325,*) '# t = ', t
        !DO k=1,lz
        !    write(325,*) k, SUM(phi(:,:,k)) 
        !END DO

        DO k=1,lz
          write(389,*) k, sum(c_plus(:,:,k))
          write(390,*) k, sum(phi(:,:,k))
        ENDDO
        ! --------------------------- Ade : 13/02/2017 ---------------------------------------------------------------
        call electrostatic_pot ! Ade: The routine is called in order to compute Phi_tot which is used in smolu
        !WRITE(387,*) '# t = ', t
        !half = lz/2
        !DO k=1,ly
        !    WRITE(387,*) k, (phi_tot(:,k,half))
        !ENDDO
        call smolu
        !WRITE(388,*) '# t = ', t
        !DO k=1,ly
        !    WRITE(388,*) k, (phi_tot(:,k,half))
        !ENDDO
        call charge_test
        write(1325,*) '# Iteration ', t
        write(1326,*) '# Iteration ', t
        write(1327,*) '# Iteration ', t
        DO k=1,lz
            write(1325,*) k, SUM(solute_force(:,:,k,1)) ! Ade : The fluid is moving in the y-direction whenever a slit 
                                                        ! case is imposed, as the walls are located at z = 0 and z = L
                                                        ! which is the reason why we are observing F_y(z). 2=>y and k=>z
            write(1326,*) k, SUM(solute_force(:,:,k,2)) 
            write(1327,*) k, SUM(solute_force(:,:,k,3)) 
        ENDDO 
        ! Ade : 19/03/17 three lines below for debugging purposes. To be removed
        do j=1,ly
            write(1328,*) j, solute_force(lx/2,j,lz/2,2)
        ENDDO

        ! Future work : Write some stuff out for postprocessing


        !--------------------------------- ADE -----------------------------------------------------------------


        !##################
        !# SINGULAR FORCE #
        !##################
        if( compensate_f_ext .and. convergence_reached_without_fext .and. t==tfext) then
            open(90,file="./output/f_ext-field_t0.dat")
            open(91,file="./output/vel-field_central_t0.dat")
            do i=1,lx
                do k=1,lz
                    WRITE(90,*) i, k, f_ext_x(i,py,k), f_ext_z(i,py,k)
                    WRITE(91,*) i, k, jx(i,py,k), jz(i,py,k)
                end do
            end do
            close(90)
            close(91)
        end if
        if( compensate_f_ext .and. convergence_reached_without_fext .and. t==tfext+1) then
            open(90,file="./output/f_ext-field_t1.dat")
            open(91,file="./output/vel-field_central_t1.dat")
            do i=1,lx
                do k=1,lz
                    WRITE(90,*) i, k, f_ext_x(i,py,k), f_ext_z(i,py,k)
                    WRITE(91,*) i, k, jx(i,py,k), jz(i,py,k)
                end do
            end do
            close(90)
            close(91)
        end if


        ! !print*,g,tock(timer(g)); g=g+1; call tick(timer(g)) !11


        !#####################
        !# check convergence #
        !#####################
        ! count the number of times the array is not zero
        !n1 = count(abs(jx_old)>1.0d-6)
        !n2 = count(abs(jy_old)>1.0d-6)
        !n3 = count(abs(jz_old)>1.0d-6)
        open(13,file="./output/l2err.dat")
        !if(n1/=0) djx = sum( abs(  (jx - jx_old)/jx_old ), mask= abs(jx_old)>1.0d-6) / real(n1,kind=dp)
        !if(n2/=0) djy = sum( abs(  (jy - jy_old)/jy_old ), mask= abs(jy_old)>1.0d-6) / real(n2,kind=dp)
        !if(n3/=0) djz = sum( abs(  (jz - jz_old)/jz_old ), mask= abs(jz_old)>1.0d-6) / real(n3,kind=dp)
        l2err = maxval([maxval(abs(jx-jx_old)), &
                        maxval(abs(jy-jy_old)), &
                        maxval(abs(jz-jz_old)) &
                       ])
        !l2err = maxval([djx,djy,djz])
        write(13,*) t, l2err

        if( l2err <= target_error .and. t>2 ) then
          convergence_reached = .true.
        else
          convergence_reached = .false.
        end if



        ! select your branch
        if(convergence_reached) then
          if( .not.convergence_reached_without_fext ) then
            convergence_reached_without_fext = .true.
          else if( convergence_reached_without_fext ) then
            convergence_reached_with_fext = .true.
          else
            print*,"ERROR: l.530 of equilibration.f90"
            print*,"=====  I did not anticipate this possibility. Review your if tree."
            stop
          end if
        end if

        !############################################
        !# Apply external contraints (f_ext) or not #
        !############################################
        if( convergence_reached ) then

          ! if you are already converged without then with f_ext then quit time loop. Stationary state is found.
          if( convergence_reached_without_fext .and. convergence_reached_with_fext .and. t>2) then
            exit ! loop over time steps

          ! if you have already converged without fext, but not yet with fext, then enable fext
          else if(convergence_reached_without_fext .and. .not.convergence_reached_with_fext) then
            tfext=t+1
            !################
            !## READ f_ext ##
            !################
            f_ext_loc = getinput%dp3("f_ext", [0._dp,0._dp,0._dp] )

            if(.not.compensate_f_ext) then ! the force is exerced everywhere with same intensity
              where(nature==fluid)
                f_ext_x = f_ext_loc(1)
                f_ext_y = f_ext_loc(2)
                f_ext_z = f_ext_loc(3)
              end where

            else if(compensate_f_ext) then ! force applied to a central particle only
                pd = getinput%int("dominika_particle_diameter",1)
                print*,"       Dominika's particle has diameter (lb units)", pd
                if( modulo(pd,2)==0 ) then
                  print*,"ERROR: l. 285 particle diameter must be odd"
                  print*,"-----  It is now",pd
                  stop
                end if

                if(modulo(lx,2)==0 .or. modulo(ly,2)==0 .or. modulo(lz,2)==0) then
                  print*,"ERROR: l.158 of equilibration.f90"
                  print*,"=====  when compensate_f_ext, there should be odd number of nodes in all directions"
                  print*,"lx, ly, lz =",lx,ly,lz
                  stop
                end if
                pdr = pd/2 ! nodes of the particle on the right (or left) of the particle center. If particle is diameter 3, we have 1 node on the left and 1 on the right, so pd=3, pdr=3/2=1

                f_ext_x = zerodp
                f_ext_y = zerodp
                f_ext_z = zerodp

                l=0 ! ADE: l counts the number of node within the particle
                err=.false.
                open(47,file="output/dominika_particle_shape.xyz")
                open(14,file="output/NodesInParticle.dat")
                ! ADE : We read the particle coordinates from lb.in
                pCoord = getinput%int3("particle_coordinates", defaultvalue=[lx/2+1,ly/2+1,lz/2+1] )
                px = pCoord(1)
                py = pCoord(2)
                pz = pCoord(3)

                do i=px-pdr,px+pdr
                  do j=py-pdr,py+pdr
                    do k=pz-pdr,pz+pdr
                      if( norm2(real([ i-(px), j-(py), k-(pz) ],dp)) > real(pd,dp)/2._dp ) cycle
                      if (nature(i,j,k)/=fluid) err=.true.
                      f_ext_x(i,j,k) = f_ext_loc(1)
                      f_ext_y(i,j,k) = f_ext_loc(2)
                      f_ext_z(i,j,k) = f_ext_loc(3)
                      l=l+1 ! ADE : One more point within the particle
                      WRITE(47,*)i,j,k ! use ListPointPlot3D[data,BoxRatios->{1,1,1}] in Mathematica to read this file
                      WRITE(14,*)l ! ADE : We write the number of lattice points occupied by the so-called particle
                    end do
                  end do
                end do
                close(14)
                close(47)
                if(err.eqv..true.) then
                  print*,"ERROR: l306 of equilibration.f90. Dominika's particle at a solid node"
                  stop
                end if

                ! ADE: We distribute the total force upon the particle evenly
                ! throughout the various nodes
                GL = getinput%int("geometryLabel", defaultvalue=0) ! if GL=-1 =>bulk case
                print*, " GL = ", GL
                ! ADE : I modified the following lines
                ! the idea is that whenever we have a slit kind of geometry,
                ! we do not want to add a compensation force in the rest of
                ! the nodes, as we believe that the force will dissipate
                ! within the walls
                if (GL==-1) then
                 where(f_ext_x==f_ext_loc(1) .and. f_ext_y==f_ext_loc(2).and.f_ext_z==f_ext_loc(3) )
                  f_ext_x = -f_ext_loc(1)/(fluid_nodes) +f_ext_x/l
                  f_ext_y = -f_ext_loc(2)/(fluid_nodes) +f_ext_y/l
                  f_ext_z = -f_ext_loc(3)/(fluid_nodes) +f_ext_z/l
                 else where
                  f_ext_x = -f_ext_loc(1)/(fluid_nodes)
                  f_ext_y = -f_ext_loc(2)/(fluid_nodes)
                  f_ext_z = -f_ext_loc(3)/(fluid_nodes)
                 end where
                   if( any(abs([sum(f_ext_x)/fluid_nodes,sum(f_ext_y)/fluid_nodes,sum(f_ext_z)/fluid_nodes])> eps ) ) then
                     print*,"ERROR: l.215 of equilibration.f90"
                     print*,"=====  The compensation is not well-implemented."
                     print*,"       sum(f_ext_x)=",sum(f_ext_x)
                     print*,"       sum(f_ext_y)=",sum(f_ext_y)
                     print*,"       sum(f_ext_z)=",sum(f_ext_z)
                     stop
                   end if
                else
                 where(f_ext_x==f_ext_loc(1) .and. f_ext_y==f_ext_loc(2) .and.f_ext_z==f_ext_loc(3) )
                  f_ext_x = f_ext_x/l
                  f_ext_y = f_ext_y/l
                  f_ext_z = f_ext_z/l
                 else where
                   f_ext_x = zerodp
                   f_ext_y = zerodp
                   f_ext_z = zerodp
                 end where
                endif
                ! ADE : end of modification

                where(nature/=fluid)
                  f_ext_x = zerodp
                  f_ext_y = zerodp
                  f_ext_z = zerodp
                end where

                print*,"       I have applied a compensating background"
            end if
          end if
        end if
        !#########################
        !# END OF SINGULAR FORCE #
        !#########################
     
        ! ADE : I added the following lines in order to write every so often the flux/velocity field values
        ! (j = rho*v), so that we can observe the trainsient time behaviour of the flow field.

        !print*, 't=',t,'mod=',mod(t,print_every) 
	if ((print_every.gt.0).and.(mod(t,print_every)==0)) then ! we divide "print_every" by the iteration time step. When remainder is zero
	                                                         ! the velocity field is written on vel-fieldTIME_*.dat (*=1,2,3,4,....)
         if( compensate_f_ext ) then
		    write(ifile,'(a,i0,a)') './output/vel-fieldTIME_', t,'.dat'
		    !print*,TRIM(ADJUSTL(ifile))
		    open(92,file=TRIM(ADJUSTL(ifile)))
		    do i=1,lx
			  do k=1,lz
			    WRITE(92,*) i, k, jx(i,py,k), jz(i,py,k)
			    !print*, i, k, jx(i,py,k), jz(i,py,k)
			  end do
		    end do
		   close(92)
         else
		    write(ifile,'(a,i0,a)') './output/vel-fieldTIME_', t,'.dat'
		    !print*,TRIM(ADJUSTL(ifile))
		    open(92,file=TRIM(ADJUSTL(ifile)))
              GL = getinput%int("geometryLabel", defaultvalue=0) ! if GL=-1 =>bulk case
              if(GL==2) then 
                do j=1,ly
                    WRITE(92,*) j, sum(jx(:,j,:)), sum(jy(:,j,:)), sum(jz(:,j,:))
                end do
              else
                do k=1,lz
                    WRITE(92,*) k, sum(jx(:,:,k)), sum(jy(:,:,k)), sum(jz(:,:,k))
                end do
              endif 
		   close(92)
		 endif
		! ----------------------------- Ade -----------------------------------------------
		! ADE : I added the following part for debugging purposes
		!write(ifile2,'(a,i0,a)') './output/c_plus_alongZTIME_', t,'.dat'
		!open(3180, file=TRIM(ADJUSTL(ifile2)))
		!DO k=supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
		!    WRITE(3180,*) k, SUM(c_plus(:,:,k))
		!ENDDO
		!close(3180)
		! ----------------------------- Ade -----------------------------------------------
    endif
        ! ADE : end of modification


  end do ! end of temporal loop


! **********************************************
! ************ Ade : POSTPROCESSING ************
!***********************************************

  ! 1. Print velocity field
  WRITE(66,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(67,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(68,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(56,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(57,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(58,*)"# Steady state with convergence criteria", REAL(target_error)

  
  GL = getinput%int("geometryLabel", defaultvalue=0) ! if GL=-1 =>bulk case
  if (GL==2) then ! Cylindrical geometry
   Jxx = 0
   Jyy = 0
   Jzz = 0
   DO i=1,lx
        Jxx = jx(i,ly/2,lz/2)
        Jyy = jy(i,ly/2,lz/2)
        Jzz = jz(i,ly/2,lz/2)
        !Jxx = jx(lx/2,j,lz/2)
        !Jyy = jy(lx/2,j,lz/2)
        !Jzz = jz(lx/2,j,lz/2)
    WRITE(67,*) i, Jxx, Jyy, Jzz
    ENDDO
        !WRITE(56,*) k, SUM(density(:,:,k),mask=node%nature==fluid)/ MAX( COUNT(density(:,:,k)>eps) ,1)
  else
    DO k=1,lz
        WRITE(66,*) k, SUM(jx(:,:,k)), SUM(jy(:,:,k)), SUM(jz(:,:,k))
        WRITE(56,*) k, SUM(density(:,:,k))/ MAX( COUNT(density(:,:,k)>eps) ,1)
    END DO
    DO k=1,ly
        WRITE(67,*) k, SUM(jx(:,k,:)), SUM(jy(:,k,:)), SUM(jz(:,k,:))
        WRITE(57,*) k, SUM(density(:,k,:))/ MAX( COUNT(density(:,k,:)>eps) ,1)
    END DO
    DO k=1,lx
        WRITE(68,*) k, SUM(jx(k,:,:)), SUM(jy(k,:,:)), SUM(jz(k,:,:))
        WRITE(58,*) k, SUM(density(k,:,:))/ MAX( COUNT(density(k,:,:)>eps) ,1)
    END DO
  endif
 
 ! 2. Solute Force
 DO k=1,lz
     write(1316,*) k, SUM(solute_force(:,:,k,1)) ! Ade : The fluid is moving in the y-direction whenever a slit 
                                                 ! case is imposed, as the walls are located at z = 0 and z = L
                                                 ! which is the reason why we are observing F_y(z). 2=>y and k=>z
     write(1323,*) k, SUM(solute_force(:,:,k,2)) 
     write(1324,*) k, SUM(solute_force(:,:,k,3)) 
 ENDDO 

 ! 3. Potential PHI
DO k=1,lz
    write(325,*) k, SUM(phi(:,:,k)) 
END DO
  close(79)
  close(80)
  CLOSE(65)
  
  CLOSE(66)
  CLOSE(67)
  CLOSE(68)
  CLOSE(56)
  CLOSE(57)
  CLOSE(58)
  CLOSE(89)
  CLOSE(90)
  CLOSE(91)
  close(316)
  close(323)
  close(324)
  close(325)
  close(1316)
  close(1323)
  close(1324)
  close(1325)
  close(1326)
  close(1327)
  close(1328)
  close(1389)
  close(1390)
  close(1391)


  !
  ! Print velocity 2D profilew
  !
  OPEN(69, FILE="output/mass-flux_field_2d_at_x.eq.1.dat")
  DO j=1,ly
      DO k=1,lz
          WRITE(69,*) j, k, jy(1,j,k), jz(1,j,k)
      END DO
  END DO
  CLOSE(69)


  if( compensate_f_ext ) then
    print*,"       The particle is located at ", px,py,pz
    open(90,file="./output/f_ext-field.dat")
    open(91,file="./output/vel-field_central.dat")
    do i=1,lx
      do k=1,lz
        WRITE(90,*) i, k, f_ext_x(i,py,k), f_ext_z(i,py,k)
        WRITE(91,*) i, k,      jx(i,py,k),      jz(i,py,k)
      end do
    end do
    close(90)
    close(91)
  end if
 close(387)
 close(388)
 close(389)

  ! put back arrays into original types
  node%solventdensity = density
  node%solventflux(x) = jx
  node%solventflux(y) = jy
  node%solventflux(z) = jz
  node%nature = nature

end subroutine equilibration
