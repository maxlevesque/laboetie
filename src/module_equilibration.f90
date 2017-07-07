module module_equilibration
    implicit none
    private
    public equilibration
contains

SUBROUTINE equilibration( jx, jy, jz)

    USE precision_kinds, only: i2b, dp, sp
    USE system, only: fluid, supercell, node, lbm, n, pbc, solute_force, c_plus, c_minus, phi, Phi_tot
    use module_collision, only: collide
    use module_input, only: getinput
    USE constants, only: x, y, z, zerodp
    USE mod_time, only: tick, tock
    USE myallocations
    use module_bounceback, only: bounceback
    use module_propagation, only: propagation
    use module_advect, only: advect

    implicit none
    real(dp), intent(inout), dimension(:,:,:) :: jx, jy, jz
    integer :: t,i,j,k,l,ip,jp,kp, lmin, lmax, timer(100), g, ng, pdr, pd, ios, px, py, pz, pCoord(3)
    integer :: fluid_nodes, print_frequency, supercellgeometrylabel, tfext, print_files_frequency, GL, print_every
    integer(kind(fluid)), allocatable, dimension(:,:,:) :: nature
    real(dp) :: n_loc, f_ext_loc(3), l2err, target_error, Jxx, Jyy, Jzz
    REAL(dp) :: vmaxx, vmaxy, vmaxz, vmax
    REAL(dp) :: vmaxx_old, vmaxy_old, vmaxz_old, vmax_old
    real(dp), allocatable, dimension(:,:,:) :: density, n_old, jx_old, jy_old, jz_old, f_ext_x, f_ext_y, f_ext_z,& 
                                               F1, F2, F3
    real(dp), allocatable, dimension(:) :: a0, a1
    integer, allocatable, dimension(:) :: cx, cy, cz
    logical :: convergence_reached, compensate_f_ext, convergence_reached_without_fext, convergence_reached_with_fext, err
    REAL(dp), PARAMETER :: eps=EPSILON(1._dp)
    LOGICAL :: write_total_mass_flux
    integer, allocatable :: il(:,:), jl(:,:), kl(:,:), l_inv(:)
    integer(i2b) :: lx, ly, lz, half
    integer(i2b) :: n1,n2,n3 ! Ade : dummy for reading purposes
    character(200) :: ifile, ifile2
    LOGICAL :: RestartPNP = .TRUE.
    integer :: maxEquilibrationTimestep

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

    RestartPNP = getinput%log("RestartPNP", .TRUE.)

    !
    ! Print info to terminal every that number of steps
    !
    print_frequency = getinput%int('print_frequency', defaultvalue=max(int(50000/(lx*ly*lz)),1), assert=">0" ) ! this number is my own optimal. To be generalized on strong criteria some day.

    !
    ! WRITE velocity profiles to terminal every that number of steps
    !
    print_files_frequency = getinput%int( "print_files_frequency", defaultvalue = huge(1) )

    print_every = getinput%int("print_every", defaultvalue=0) ! reads from lb.in file
                                                                 ! the frequency of printing time
                                                                 ! the default value needs to be changed eventually
    fluid_nodes = count( node%nature==fluid )

    ! Max : I had 1.D-8 before ADE's modification (June 21)
    target_error = getinput%dp("target_error", 1.D-10)

    allocate( density(lx,ly,lz), source=sum(n,4), stat=ios)
    if (ios /= 0) stop "density: Allocation request denied"

    allocate( l_inv(lmin:lmax) , stat=ios)
    if (ios /= 0) stop "l_inv: Allocation request denied"
    do l = lmin, lmax
        l_inv(l) = lbm%vel(l)%inv
    end do

    if(.not.allocated(solute_force)) allocate(solute_force(lx,ly,lz,x:z),source=0.0_dp)
    allocate( jx_old (lx,ly,lz), source=0._dp )
    allocate( jy_old (lx,ly,lz), source=0._dp )
    allocate( jz_old (lx,ly,lz), source=0._dp )

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
    
    compensate_f_ext = getinput%log( "compensate_f_ext", defaultvalue = .false.)
    if(compensate_f_ext) then
        open(79, file = "./output/v_centralnode.dat")
        open(80, file = "./output/rho_centralnode.dat")
    endif

    write_total_mass_flux = getinput%log( "write_total_mass_flux", defaultvalue = .false.)
        if( write_total_mass_flux ) open(65, file = "./output/total_mass_flux.dat" )



    PRINT*
    PRINT*,'Lattice Boltzmann'
    PRINT*,'================='
    PRINT*,'       step'
    PRINT*,'       ----'


    ! ADE : We initialise tfext, which is the time from when the force f_ext is applied
    ! upon a certain number of nodes. 
    tfext = HUGE(tfext)
    ! We also init l2err, the convergence error
    l2err = -1


    maxEquilibrationTimestep = getinput%int( 't_equil' , defaultvalue = -1)
    ! We start by equilibrating the densities, fluxes and other moments without the external forces.
    ! This equilibration step is done up to equilibration 
    ! or if the timestep t gets higher than a maximum value called maxEquilibrationTimestep


    !
    ! TIME STEPS (in lb units)
    !
    do t = 1, huge(t)

        if( t < maxEquilibrationTimestep ) then
            f_ext_x = zerodp
            f_ext_y = zerodp
            f_ext_z = zerodp
        else
            f_ext_x = f_ext_loc(1) ! Max: as I understand this, f_ext_loc(1:3)=0 but if convergence is reached. Once convergence is reached, we read f_ext_loc from input file.
            f_ext_y = f_ext_loc(2)
            f_ext_z = f_ext_loc(3)
        endif

        !
        ! Print sdtout timestep, etc.
        !
        if( modulo(t, print_frequency) == 0) PRINT*, t, real(l2err),"(target",real(target_error,4),")"

        !
        ! Print velocity profiles, etc.
        !
        if( modulo(t, print_files_frequency) == 0 ) then
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
        end if

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


        !###############
        !# Bounce Back # to simplify propagation
        !###############
        call bounceback(n, nature)


        !###############
        !# PROPAGATION #
        !###############
        call propagation(n, lmin, lmax, lx, ly, lz, il, jl, kl)

        !###############
        !# CHECK POPULATIONS ARE POSITIVE
        !###############
        IF( ANY(n<0) ) ERROR STOP "In equilibration, the population n(x,y,z,vel) < 0"

        !
        ! Finally update densities after the propagation and check it
        ! Densities are the sum of all velocities of a local population
        !
        density = SUM(n,4)

        ! WRITE the total density
        IF( write_total_mass_flux ) WRITE(65,*) t, REAL([  SUM(jx), SUM(jy), SUM(jz)  ])

        !
        ! backup moment density (velocities) to test convergence at the end of the timestep
        !
        jx_old = jx
        jy_old = jy
        jz_old = jz

        call update_solventCurrent( jx, jy, jz, n, cx, cy, cz, F1, F2, F3)

        call advect( density, jx, jy, jz )

        call sor    ! compute the electric potential phi with the Successive OverRelation method (SOR)

        call electrostatic_pot ! Ade: The routine is called in order to compute Phi_tot which is used in smolu

        call smolu

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
            write(1328,*) j, solute_force(max(lx/2,1),j,max(lz/2,1),2)
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


        !#####################
        !# check convergence #
        !#####################
        open(13,file="./output/l2err.dat")
        l2err = maxval([maxval(abs(jx-jx_old)), &
                        maxval(abs(jy-jy_old)), &
                        maxval(abs(jz-jz_old)) &
                       ])
        write(13,*) t, l2err

! convergence_reached = .false. ! by default
! if( t > 2 ) then
!     if( check_convergence( jx, jx_old ) ) then
!         if( check_convergence( jy, jy_old ) ) then
!             if( check_convergence( jz, jz_old ) ) then
!                 convergence_reached = .true.
!             end if
!         end if
!     end if
! end if

! pure function check_convergence( A, A_old, criteria )
!     implicit none
!     logical :: convergence_reached
!     real(dp) :: l2err
!     if( maxval(abs(A-A_old)) <= criteria ) check_convergence = .true.
! end function check_convergence

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
            f_ext_loc = getinput%dp3("f_ext", defaultvalue= [0._dp,0._dp,0._dp] )

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
            end if 
            close(92)
        end if
    endif

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
        Jxx = jx(i,max(ly/2,1),max(lz/2,1))
        Jyy = jy(i,max(ly/2,1),max(lz/2,1))
        Jzz = jz(i,max(ly/2,1),max(lz/2,1))
        WRITE(67,*) i, Jxx, Jyy, Jzz
    END DO
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
end if
  CLOSE(68)
  CLOSE(58)
    CLOSE(67)
    CLOSE(57)
    CLOSE(66)
    CLOSE(56)
 
 ! 2. Solute Force
 DO k=1,lz
     write(1316,*) k, SUM(solute_force(:,:,k,1)) ! Ade : The fluid is moving in the y-direction whenever a slit 
                                                 ! case is imposed, as the walls are located at z = 0 and z = L
                                                 ! which is the reason why we are observing F_y(z). 2=>y and k=>z
     write(1323,*) k, SUM(solute_force(:,:,k,2)) 
     write(1324,*) k, SUM(solute_force(:,:,k,3)) 
 ENDDO 
  close(1316)
  close(1323)
  close(1324)

 ! 3. Potential PHI
DO k=1,lz
    write(325,*) k, SUM(phi(:,:,k)) 
END DO
  close(325)

  close(79)
  close(80)
  CLOSE(65)
  CLOSE(89)
  CLOSE(90)
  CLOSE(91)
  close(316)
  close(323)
  close(324)
 close(387)
 close(388)
 close(389)
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

end subroutine equilibration

subroutine update_solventCurrent( jx, jy, jz, n, cx, cy, cz, F1, F2, F3)
    use precision_kinds, only: dp
    implicit none
    real(dp), intent(inout), dimension(:,:,:) :: jx, jy, jz
    real(dp), intent(in) :: n(:,:,:,:), F1(:,:,:), F2(:,:,:), F3(:,:,:)
    integer, intent(in) :: cx(:), cy(:), cz(:)
    integer :: lmin, lmax, l
    lmin = lbound(cx,1)
    lmax = ubound(cx,1)
    ! update momentum densities after the propagation
    ! this is completely local in space and my be parallelized very well
    jx = F1/2._dp
    jy = F2/2._dp
    jz = F3/2._dp
    !$OMP PARALLEL DO DEFAULT(NONE)&
    !$OMP PRIVATE(l)&
    !$OMP SHARED(lmin,lmax,n,cx,cy,cz)&
    !$OMP REDUCTION(+:jx)&
    !$OMP REDUCTION(+:jy)&
    !$OMP REDUCTION(+:jz)
    do l=lmin,lmax
        jx = jx +n(:,:,:,l)*cx(l)
        jy = jy +n(:,:,:,l)*cy(l)
        jz = jz +n(:,:,:,l)*cz(l)
    end do
    !$OMP END PARALLEL DO
end subroutine update_solventCurrent

end module module_equilibration
