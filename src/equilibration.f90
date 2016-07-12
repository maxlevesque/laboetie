SUBROUTINE equilibration

    USE precision_kinds, only: i2b, dp, sp
    USE system, only: fluid, supercell, node, lbm, n, pbc
    use module_collision, only: collide
    use module_input, only: getinput
    USE constants, only: x, y, z, zerodp
    USE mod_time, only: tick, tock

    implicit none
    integer :: t,i,j,k,l,ip,jp,kp,n1,n2,n3, lmin, lmax, timer(100), g, ng, pdr, pd, ios, px, py, pz, pCoord(3)
    integer :: fluid_nodes, print_frequency, supercellgeometrylabel, tfext, print_files_frequency, GL
    integer(kind(fluid)), allocatable, dimension(:,:,:) :: nature
    real(dp) :: n_loc, f_ext_loc(3), l2err, target_error
    REAL(dp) :: vmaxx, vmaxy, vmaxz, vmax
    REAL(dp) :: vmaxx_old, vmaxy_old, vmaxz_old, vmax_old
    real(dp), allocatable, dimension(:,:,:) :: density, jx, jy, jz, n_old, jx_old, jy_old, jz_old, f_ext_x, f_ext_y, f_ext_z
    real(dp), allocatable, dimension(:) :: a0, a1
    integer, allocatable, dimension(:) :: cx, cy, cz
    logical :: convergence_reached, compensate_f_ext, convergence_reached_without_fext, convergence_reached_with_fext, err
    REAL(dp), PARAMETER :: eps=EPSILON(1._dp)
    LOGICAL :: write_total_mass_flux
    integer, allocatable :: il(:,:), jl(:,:), kl(:,:), l_inv(:)

    !
    ! laboetie doesnt work for charged solutes
    !
    IF( ABS(getinput%dp('sigma', zerodp)) > eps ) THEN
        print*,"ERROR: laboetie can only consider uncharged systems."
        print*,"===== Dont tell Benjamin you'd like to see such feature in Laboetie :)"
        print*,"Hi Benjamin. I'm sure it is you testing this! grrrr :))"
        ERROR STOP
    END IF

    lmin = lbm%lmin
    lmax = lbm%lmax

    supercellgeometrylabel = supercell%geometry%label ! -1 for solid free cell

    n1 = getinput%int("lx", assert=">0")
    n2 = getinput%int("ly", assert=">0")
    n3 = getinput%int("lz", assert=">0")

    !
    ! Print info to terminal every that number of steps
    !
    print_frequency = getinput%int('print_frequency', defaultvalue=max(int(50000/(n1*n2*n3)),1), assert=">0" ) ! this number is my own optimal. To be generalized on strong criteria some day.

    !
    ! WRITE velocity profiles to terminal every that number of steps
    !
    print_files_frequency = getinput%int("print_files_frequency", HUGE(1))

    fluid_nodes = count( node%nature==fluid )

    ! Max : I had 1.D-8 before ADE's modification (June 21)
    target_error = getinput%dp("target_error", 1.D-10)

    allocate( density(n1,n2,n3), source=node%solventdensity, stat=ios)
    if (ios /= 0) stop "density: Allocation request denied"

    allocate( l_inv(lmin:lmax) , stat=ios)
    if (ios /= 0) stop "l_inv: Allocation request denied"
    do l = lmin, lmax
        l_inv(l) = lbm%vel(l)%inv
    end do


    allocate( jx     (n1,n2,n3), source=node%solventflux(x))
    allocate( jx_old (n1,n2,n3) )
    allocate( jy     (n1,n2,n3), source=node%solventflux(y))
    allocate( jy_old (n1,n2,n3) )
    allocate( jz     (n1,n2,n3), source=node%solventflux(z))
    allocate( jz_old (n1,n2,n3) )
    jx = 0
    jy = 0
    jz = 0
    jx_old = 0
    jy_old = 0
    jz_old = 0

    OPEN(66, FILE="output/mass-flux_profile_along_z.dat")
    OPEN(67, FILE="output/mass-flux_profile_along_y.dat")
    OPEN(68, FILE="output/mass-flux_profile_along_x.dat")
    WRITE(66,*) "# z, <ρ.v_x>_{x,y}, <ρ.v_y>_{x,y}, <ρ.v_z>_{x,y}"
    WRITE(67,*) "# y, <ρ.v_x>_{x,z}, <ρ.v_y>_{x,z}, <ρ.v_z>_{x,z}"
    WRITE(68,*) "# x, <ρ.v_x>_{y,z}, <ρ.v_y>_{y,z}, <ρ.v_z>_{y,z}"

    OPEN(56, FILE="output/mean-density_profile_along_z.dat")
    OPEN(57, FILE="output/mean-density_profile_along_y.dat")
    OPEN(58, FILE="output/mean-density_profile_along_x.dat")

    allocate( nature (n1,n2,n3), source=node%nature)
    allocate( f_ext_x(n1,n2,n3), source=zerodp)
    allocate( f_ext_y(n1,n2,n3), source=zerodp)
    allocate( f_ext_z(n1,n2,n3), source=zerodp)

    f_ext_loc = zerodp ! this is important and spagetty like... please read carefuln2 before modifying this line
    allocate( cx(lmax), source=lbm%vel(:)%coo(1))
    allocate( cy(lmax), source=lbm%vel(:)%coo(2))
    allocate( cz(lmax), source=lbm%vel(:)%coo(3))
    allocate( a0(lmax), source=lbm%vel(:)%a0)
    allocate( a1(lmax), source=lbm%vel(:)%a1)

    !
    ! Tabulate the index of the node one finishes if one starts from a node and a velocity index l
    ! per direction
    !
    allocate( il(lbm%lmin:lbm%lmax, 1:n1), stat=ios)
    if (ios /= 0) stop "il: Allocation request denied"
    allocate( jl(lbm%lmin:lbm%lmax, 1:n2), stat=ios)
    if (ios /= 0) stop "jl: Allocation request denied"
    allocate( kl(lbm%lmin:lbm%lmax, 1:n3), stat=ios)
    if (ios /= 0) stop "kl: Allocation request denied"
    do l= lmin, lmax
        il(l,:) = [( pbc(i+cx(l),x) ,i=1,n1 )]
        jl(l,:) = [( pbc(j+cy(l),y) ,j=1,n2 )]
        kl(l,:) = [( pbc(k+cz(l),z) ,k=1,n3 )]
    end do

    convergence_reached_without_fext = .false.
    convergence_reached_with_fext = .false.
    compensate_f_ext = getinput%log("compensate_f_ext",.false.)
    if(compensate_f_ext) open(79,file="./output/v_centralnode.dat")

    write_total_mass_flux = getinput%log("write_total_mass_flux", .FALSE.)
    IF( write_total_mass_flux ) THEN
        OPEN( 65, FILE="output/total_mass_flux.dat" )
    END IF



    PRINT*
    PRINT*,'Lattice Boltzmann'
    PRINT*,'================='
    PRINT*,'       step     flux max           error         target error'
    PRINT*,'       ----------------------------------------------------------------------------------'


    !
    ! TIME STEPS (in lb units)
    !
    do t=1,HUGE(t)

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
            DO k=1,n3
                WRITE(56,*) k, SUM(density(:,:,k))/ MAX( COUNT(density(:,:,k)>eps)  ,1)
                WRITE(66,*) k, SUM(jx(:,:,k)), SUM(jy(:,:,k)), SUM(jz(:,:,k))
            END DO
            DO k=1,n2
                WRITE(57,*) k, SUM(density(:,k,:))/ MAX( COUNT(density(:,k,:)>eps)  ,1)
                WRITE(67,*) k, SUM(jx(:,k,:)), SUM(jy(:,k,:)), SUM(jz(:,k,:))
            END DO
            DO k=1,n1
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
        end if

        !
        ! backup moment density (velocities) to test convergence at the end of the timestep
        !
        jx_old = jx
        jy_old = jy
        jz_old = jz

        !print*,g,tock(timer(g)); g=g+1; call tick(timer(g)) !3

        !
        ! Collision step
        !
        call collide(n, density, jx, jy, jz, f_ext_x, f_ext_y, f_ext_z)

        ! print velocity profile if you need/want it
        ! if( modulo(t, print_frequency) == 0) then
        !    call velocity_profiles(t) ! print velocity profiles
        ! end if

        !
        ! Bounce back (boundpm) to simplify propagation step
        !
        do concurrent(l=lmin:lmax:2)
            do concurrent(k=1:n3)
                kp = kl(l,k)
                !kp=pbc(k+cz(l),z)
                do concurrent(j=1:n2)
                    jp = jl(l,j)
                    !jp=pbc(j+cy(l),y)
                    do concurrent(i=1:n1)
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

        !
        ! propagation
        !
        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(n,n3,n2,n1,lmin,lmax,cz,cy,cx,il,jl,kl) &
        !$OMP PRIVATE(l,k,j,i,ip,jp,kp,n_old)
        do l=lmin,lmax
            n_old = n(:,:,:,l)
            do k=1,n3
                kp = kl(l,k)
                do j=1,n2
                    jp = jl(l,j)
                    do i=1,n1
                        ip = il(l,i)
                        n(ip,jp,kp,l) = n_old(i,j,k)
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        !
        ! The populations, since they are probabilities, must never be negative
        !
        IF( ANY(n<0) ) ERROR STOP "In equilibration_new, the population n(x,y,z,vel) < 0"

        !
        ! Update densities after the propagation and check it
        ! Densities are the sum of all velocities of a local population
        !
        density = SUM(n,4)

        !
        ! WRITE the total density
        !
        IF( write_total_mass_flux ) THEN
            WRITE(65,*) t, REAL([  SUM(jx), SUM(jy), SUM(jz)  ])
        END IF


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
        jx=0
        jy=0
        jz=0
        do l=lmin,lmax
            jx = jx +n(:,:,:,l)*cx(l)
            jy = jy +n(:,:,:,l)*cy(l)
            jz = jz +n(:,:,:,l)*cz(l)
        end do

        !
        ! Dominika
        !
        if( compensate_f_ext .and. convergence_reached_without_fext .and. t==tfext) then
            open(90,file="./output/f_ext-field_t0.dat")
            open(91,file="./output/vel-field_central_t0.dat")
            do i=1,n1
                do k=1,n3
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
            do i=1,n1
                do k=1,n3
                    WRITE(90,*) i, k, f_ext_x(i,py,k), f_ext_z(i,py,k)
                    WRITE(91,*) i, k, jx(i,py,k), jz(i,py,k)
                end do
            end do
            close(90)
            close(91)
        end if


        ! !print*,g,tock(timer(g)); g=g+1; call tick(timer(g)) !11


        !
        ! check convergence
        !
        open(13,file="./output/l2err.dat")
        l2err = maxval([maxval(abs(jx-jx_old)), &
                        maxval(abs(jy-jy_old)), &
                        maxval(abs(jz-jz_old)) &
                       ])
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
            print*,"ERROR: l.376 of equilibration_new.f90"
            print*,"=====  I did not anticipate this possibility. Review your if tree."
            stop
          end if
        end if

        !
        ! Apply external contraints (f_ext) or not
        !
        if( convergence_reached ) then

          ! if you are already converged without then with f_ext then quit time loop. Stationary state is found.
          if( convergence_reached_without_fext .and. convergence_reached_with_fext .and. t>2) then
            exit ! loop over time steps

          ! if you have already converged without fext, but not yet with fext, then enable fext
          else if(convergence_reached_without_fext .and. .not.convergence_reached_with_fext) then
            tfext=t+1
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

                if(modulo(n1,2)==0 .or. modulo(n2,2)==0 .or. modulo(n3,2)==0) then
                  print*,"ERROR: l.158 of equilibration_new.f90"
                  print*,"=====  when compensate_f_ext, there should be odd number of nodes in all directions"
                  print*,"n1, n2, n3 =",n1,n2,n3
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
                pCoord = getinput%int3("particle_coordinates", defaultvalue=[n1/2+1,n2/2+1,n3/2+1] )
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
                  print*,"ERROR: l306 of equilibration_new.f90. Dominika's particle at a solid node"
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
                     print*,"ERROR: l.215 of equilibration_new.f90"
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

  end do

  close(79)
  CLOSE(65)

  !
  ! Print velocity 1D velocity field
  !
  WRITE(66,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(67,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(68,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(56,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(57,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(58,*)"# Steady state with convergence criteria", REAL(target_error)
  DO k=1,n3
      WRITE(66,*) k, SUM(jx(:,:,k)), SUM(jy(:,:,k)), SUM(jz(:,:,k))
      WRITE(56,*) k, SUM(density(:,:,k))/ MAX( COUNT(density(:,:,k)>eps) ,1)
  END DO
  DO k=1,n2
      WRITE(67,*) k, SUM(jx(:,k,:)), SUM(jy(:,k,:)), SUM(jz(:,k,:))
      WRITE(57,*) k, SUM(density(:,k,:))/ MAX( COUNT(density(:,k,:)>eps) ,1)
  END DO
  DO k=1,n1
      WRITE(68,*) k, SUM(jx(k,:,:)), SUM(jy(k,:,:)), SUM(jz(k,:,:))
      WRITE(58,*) k, SUM(density(k,:,:))/ MAX( COUNT(density(k,:,:)>eps) ,1)
  END DO
  CLOSE(66)
  CLOSE(67)
  CLOSE(68)
  CLOSE(56)
  CLOSE(57)
  CLOSE(58)

  !
  ! Print velocity 2D profilew
  !
  OPEN(69, FILE="output/mass-flux_field_2d_at_x.eq.1.dat")
  DO j=1,n2
      DO k=1,n3
          WRITE(69,*) j, k, jy(1,j,k), jz(1,j,k)
      END DO
  END DO
  CLOSE(69)


  if( compensate_f_ext ) then
    print*,"       The particle is located at ", px,py,pz
    open(90,file="./output/f_ext-field.dat")
    open(91,file="./output/vel-field_central.dat")
    do i=1,n1
      do k=1,n3
        WRITE(90,*) i, k, f_ext_x(i,py,k), f_ext_z(i,py,k)
        WRITE(91,*) i, k,      jx(i,py,k),      jz(i,py,k)
      end do
    end do
    close(90)
    close(91)
  end if

  ! put back arrays into original types
  node%solventdensity = density
  node%solventflux(x) = jx
  node%solventflux(y) = jy
  node%solventflux(z) = jz
  node%nature = nature

end subroutine equilibration
