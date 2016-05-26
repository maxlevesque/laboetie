SUBROUTINE equilibration

    USE precision_kinds, only: i2b, dp, sp
    USE system, only: fluid, supercell, node, lbm, n, pbc
    USE populations, only: update_populations
    use module_input, only: getinput
    USE constants, only: x, y, z, zerodp
    USE mod_time, only: tick, tock

    implicit none
    integer :: t,i,j,k,l,ip,jp,kp,n1,n2,n3, lmin, lmax, timer(100), g, ng, pdr, pd, ios
    integer :: fluid_nodes, print_frequency, supercellgeometrylabel, tfext, print_files_frequency
    integer(kind(fluid)), allocatable, dimension(:,:,:) :: nature
    real(dp) :: n_loc, f_ext_loc(3), l2err, target_error
    real(dp) :: vmaxx, vmaxy, vmaxz, vmax
    real(dp) :: vmaxx_old, vmaxy_old, vmaxz_old, vmax_old
    real(dp), allocatable, dimension(:,:,:) :: density, jx, jy, jz, n_old, jx_old, jy_old, jz_old, f_ext_x, f_ext_y, f_ext_z
    real(dp), allocatable, dimension(:) :: a0, a1
    integer, allocatable, dimension(:) :: cx, cy, cz
    logical :: convergence_reached, compensate_f_ext, convergence_reached_without_fext, convergence_reached_with_fext, err
    real(dp), PARAMETER :: eps=EPSILON(1._dp)
    logical :: write_total_mass_flux
    integer, allocatable :: il(:,:), jl(:,:), kl(:,:), l_inv(:)
    integer :: px, py, pz, pCoord(3)

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

    target_error = getinput%dp("target_error", 1.D-8)

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
        IF( MODULO(t, print_frequency) == 0) THEN
            vmaxx = MAXVAL( ABS(jx)/density, density>eps)
            vmaxy = MAXVAL( ABS(jy)/density, density>eps)
            vmaxz = MAXVAL( ABS(jz)/density, density>eps)
            vmax  = MAX( vmaxx, vmaxy, vmaxz )
            PRINT*, t, real( [vmax, l2err, target_error] ,sp)
        END IF

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
          ! We read the particle coordinates from lb.in
          pCoord = getinput%int3( "particle_coordinates", defaultvalue=[n1/2+1,n2/2+1,n3/2+1] )
          px = pCoord(1)
          py = pCoord(2)
          pz = pCoord(3)
          write(79,*) t-tfext, jx(px,py,pz), jy(px,py,pz), jz(px,py,pz)
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
        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(n,lmin,lmax,cx,jx,f_ext_x,cy,jy,f_ext_y,cz,jz,f_ext_z,a0,density,a1)&
        !$OMP PRIVATE(l)
        do l=lmin,lmax
            n(:,:,:,l) = a0(l)*density(:,:,:) + a1(l)*(cx(l)*(jx+f_ext_x)+cy(l)*(jy+f_ext_y)+cz(l)*(jz+f_ext_z))
        end do
        !$OMP END PARALLEL DO

        ! do concurrent(i=1:n1, j=1:n2, k=1:n3, l=lmin:lmax)
        !   n(i,j,k,l) = a0(l)*density(i,j,k) +a1(l)*(&
        !     cx(l)*(jx(i,j,k)+f_ext_x(i,j,k)) + cy(l)*(jy(i,j,k)+f_ext_y(i,j,k)) + cz(l)*(jz(i,j,k)+f_ext_z(i,j,k)))
        ! end do

        ! print velocity profile if you need/want it
        ! if( modulo(t, print_frequency) == 0) then
        !    call velocity_profiles(t) ! print velocity profiles
        ! end if

        !print*,g,tock(timer(g)); g=g+1; call tick(timer(g)) !5

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


        ! print*,sum(density)
        ! if( abs(sum(density)/real(n1*n2*n3,kind(density)) -1._dp) > epsilon(1._dp) ) then
        !   stop "otto"
        ! end if

        !print*,g,tock(timer(g)); g=g+1; call tick(timer(g)) !9

        ! update momentum densities after the propagation
        ! this is completely local in space and my be parallelized very well
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
        jx=jx/2
        jy=jy/2
        jz=jz/2
        ! do concurrent (i=1:n1, j=1:n2, k=1:n3)
        !   jx(i,j,k) = (jx(i,j,k) + sum(n(i,j,k,:)*cx(:)))/2._dp
        !   jy(i,j,k) = (jy(i,j,k) + sum(n(i,j,k,:)*cy(:)))/2._dp
        !   jz(i,j,k) = (jz(i,j,k) + sum(n(i,j,k,:)*cz(:)))/2._dp
        ! end do

        !
        ! Dominika
        !
        if( compensate_f_ext .and. convergence_reached_without_fext .and. t==tfext) then
            open(90,file="./output/f_ext-field_t0.dat")
            open(91,file="./output/vel-field_central_t0.dat")
            do i=1,n1
                do k=1,n3
                    WRITE(90,*) i, k, f_ext_x(i,n2/2+1,k), f_ext_z(i,n2/2+1,k)
                    WRITE(91,*) i, k,      jx(i,n2/2+1,k)     ,      jz(i,n2/2+1,k)
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
                    WRITE(90,*) i, k, f_ext_x(i,n2/2+1,k), f_ext_z(i,n2/2+1,k)
                    WRITE(91,*) i, k, jx(i,n2/2+1,k), jz(i,n2/2+1,k)
                end do
            end do
            close(90)
            close(91)
        end if


        ! !print*,g,tock(timer(g)); g=g+1; call tick(timer(g)) !11


        !
        ! check convergence
        !
        l2err = norm2(jx-jx_old+jy-jy_old+jz-jz_old)
        l2err = sqrt(  norm2(jx-jx_old)**2  +  norm2(jy-jy_old)**2  +  norm2(jz-jz_old)**2  )/fluid_nodes
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
            print*,"ERROR: l.182 of equilibration_new.f90"
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

                l=0
                err=.false.
                open(47,file="output/dominika_particle_shape.xyz")
                open(14,file="output/NodesInParticle.dat")
                do i= px-pdr, px+pdr
                  do j= py-pdr, py+pdr
                    do k= pz-pdr, pz+pdr
                      if( norm2(real([ i-px, j-py, k-pz ],dp)) > real(pd,dp)/2._dp ) cycle
                      if (nature(i,j,k)/=fluid) err=.true.
                      f_ext_x(i,j,k) = f_ext_loc(1)
                      f_ext_y(i,j,k) = f_ext_loc(2)
                      f_ext_z(i,j,k) = f_ext_loc(3)
                      l=l+1 ! One more point within the particle
                      write(47,*) i,j,k ! use ListPointPlot3D[data,BoxRatios->{1,1,1}]
                                        ! in Mathematica to read this file
                      write(14,*) l ! We write the number of lattice points occupied by the so-called particle
                    end do
                  end do
                end do
                close(47)
                if(err.eqv..true.) then
                  print*,"ERROR: l306 of equilibration_new.f90."
                  print*,"=====  The particle is on top of a solid node"
                  stop
                end if

                !
                ! We distribute the total force upon the particle evenly throughout the various nodes.
                ! The so-called particle is sensitive to the total number of nodes in the physical system.
                !
                where(f_ext_x==f_ext_loc(1) .and. f_ext_y==f_ext_loc(2) .and. f_ext_z==f_ext_loc(3) )
                  f_ext_x = -f_ext_loc(1)/(fluid_nodes) +f_ext_x/l
                  f_ext_y = -f_ext_loc(2)/(fluid_nodes) +f_ext_y/l
                  f_ext_z = -f_ext_loc(3)/(fluid_nodes) +f_ext_z/l
                else where
                  f_ext_x = -f_ext_loc(1)/(fluid_nodes)
                  f_ext_y = -f_ext_loc(2)/(fluid_nodes)
                  f_ext_z = -f_ext_loc(3)/(fluid_nodes)
                end where

                !
                ! Be sure you don't apply a force on a solid node
                !
                where(nature/=fluid)
                  f_ext_x = zerodp
                  f_ext_y = zerodp
                  f_ext_z = zerodp
                end where

                ! check that we have a compensating background, i.e., that total force is zero
                if( any( abs([sum(f_ext_x),sum(f_ext_y),sum(f_ext_z)]) > eps ) ) then
                  print*,"ERROR: l.215 of equilibration_new.f90"
                  print*,"=====  The compensation is not well-implemented."
                  print*,"       sum(f_ext_x)=",sum(f_ext_x)
                  print*,"       sum(f_ext_y)=",sum(f_ext_y)
                  print*,"       sum(f_ext_z)=",sum(f_ext_z)
                  stop
                end if

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
    print*,"       Central node is at",n1/2+1,n2/2+1,n3/2+1
    open(90,file="./output/f_ext-field.dat")
    open(91,file="./output/vel-field_central.dat")
    do i=1,n1
      do k=1,n3
        WRITE(90,*) i, k, f_ext_x(i,n2/2+1,k), f_ext_z(i,n2/2+1,k)
        WRITE(91,*) i, k,      jx(i,n2/2+1,k),      jz(i,n2/2+1,k)
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
