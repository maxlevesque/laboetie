SUBROUTINE equilibration

    USE precision_kinds, only: i2b, dp, sp
    USE system, only: fluid, supercell, node, lbm, n, pbc
    USE populations, only: update_populations
    USE input, only: input_dp3, input_dp, input_int, input_log
    USE constants, only: x, y, z, zerodp
    USE mod_time, only: tick, tock

    implicit none
    integer :: t,i,j,k,l,ip,jp,kp,n1,n2,n3,lmax,l_inv,timer(100),g,ng,pdr,pd
    integer :: fluid_nodes, print_frequency, supercellgeometrylabel, tfext, print_files_frequency
    integer(kind(fluid)), allocatable, dimension(:,:,:) :: nature
    real(dp) :: n_loc, f_ext_loc(3), l2err, target_error
    REAL(dp) :: vmaxx, vmaxy, vmaxz, vmax
    REAL(dp) :: vmaxx_old, vmaxy_old, vmaxz_old, vmax_old
    real(dp), allocatable, dimension(:,:,:) :: density, jx, jy, jz, old_n, jx_old, jy_old, jz_old, f_ext_x, f_ext_y, f_ext_z
    real(dp), allocatable, dimension(:) :: a0, a1
    integer, allocatable, dimension(:) :: cx, cy, cz
    logical :: convergence_reached, compensate_f_ext, convergence_reached_without_fext, convergence_reached_with_fext, err
    INTEGER :: lx, ly, lz
    REAL(dp), PARAMETER :: eps=EPSILON(1._dp)
    LOGICAL :: write_total_mass_flux

    !
    ! laboetie doesnt work for charged solutes
    !
    IF( ABS(input_dp('sigma', zerodp)) > eps ) THEN
        print*,"ERROR: laboetie can only consider uncharged systems."
        print*,"===== Dont tell Benjamin you'd like to see such feature in Laboetie :)"
        print*,"Hi Benjamin. I'm sure it is you testing this! grrrr :))"
        ERROR STOP
    END IF

    supercellgeometrylabel = supercell%geometry%label ! -1 for solid free cell
    lx = input_int( "lx" )
    ly = input_int( "ly" )
    lz = input_int( "lz" )

    !
    ! Print info to terminal every that number of steps 
    !
    print_frequency = input_int('print_frequency', INT(50000/(lx*ly*lz)) )
    IF( print_frequency==0) print_frequency=1

    !
    ! WRITE velocity profiles to terminal every that number of steps
    !
    print_files_frequency = input_int("print_files_frequency", HUGE(1))

    fluid_nodes = count( node%nature==fluid )

    target_error = input_dp("target_error", 1.D-10)
     
    n1 = supercell%geometry%dimensions%indicemax(1)
    n2 = supercell%geometry%dimensions%indicemax(2)
    n3 = supercell%geometry%dimensions%indicemax(3)
    allocate( density(n1,n2,n3), source=node%solventdensity )

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

    f_ext_loc = zerodp ! this is important and spagetty like... please read carefully before modifying this line
    If( lbm%lmin /= 1) error stop "lbm%lmin expected to be 1"
    lmax = lbm%lmax
    allocate( cx(lmax), source=lbm%vel(:)%coo(1))
    allocate( cy(lmax), source=lbm%vel(:)%coo(2))
    allocate( cz(lmax), source=lbm%vel(:)%coo(3))
    allocate( a0(lmax), source=lbm%vel(:)%a0)
    allocate( a1(lmax), source=lbm%vel(:)%a1)

    convergence_reached_without_fext = .false.
    convergence_reached_with_fext = .false.
    compensate_f_ext = input_log("compensate_f_ext",.false.)
    if(compensate_f_ext) open(79,file="./output/v_centralnode.dat")

    write_total_mass_flux = input_log("write_total_mass_flux", .FALSE.)
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
        ! Dominika
        !
        if( compensate_f_ext .and. convergence_reached_without_fext) then
          WRITE(79,*)t-tfext, jz(n1/2+1,n2/2+1,n3/2+1)
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
        !$OMP SHARED(n,lmax,cx,jx,f_ext_x,cy,jy,f_ext_y,cz,jz,f_ext_z,a0,density,a1)&
        !$OMP PRIVATE(l)
        do l=1,lmax
            n(:,:,:,l) = a0(l)*density &
                + a1(l)*(cx(l)*(jx+f_ext_x)+cy(l)*(jy+f_ext_y)+cz(l)*(jz+f_ext_z))
        end do
        !$OMP END PARALLEL DO

        ! do concurrent(i=1:n1, j=1:n2, k=1:n3, l=1:lmax)
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
        do concurrent(l=1:lmax:2)
            l_inv = lbm%vel(l)%inv
            do concurrent(k=1:n3)
                kp=pbc(k+cz(l),z)
                do concurrent(j=1:n2)
                    jp=pbc(j+cy(l),y)
                    do concurrent(i=1:n1)
                       ip=pbc(i+cx(l),x)
                       if( nature(i,j,k) /= nature(ip,jp,kp) ) then
                           n_loc = n(i,j,k,l)
                           n(i,j,k,l) = n(ip,jp,kp,l_inv)
                           n(ip,jp,kp,l_inv) = n_loc
                       end if
                    end do
                end do
            end do
        end do

        !
        ! propagation
        !
        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(n,n3,n2,n1,lmax,cz,cy,cx) &
        !$OMP PRIVATE(l,k,j,i,ip,jp,kp,old_n)
        do l=1,lmax
          old_n = n(:,:,:,l)
          do k=1,n3
            kp=pbc(k+cz(l),z)
            do j=1,n2
              jp=pbc(j+cy(l),y)
              do i=1,n1
                ip=pbc(i+cx(l),x)
                n(ip,jp,kp,l) = old_n(i,j,k)
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
        !$OMP SHARED(lmax,n,cx,cy,cz)&
        !$OMP REDUCTION(+:jx)&
        !$OMP REDUCTION(+:jy)&
        !$OMP REDUCTION(+:jz)
        do l=1,lmax
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
            f_ext_loc = input_dp3("f_ext", [0._dp,0._dp,0._dp] )

            if(.not.compensate_f_ext) then ! the force is exerced everywhere with same intensity
              where(nature==fluid)
                f_ext_x = f_ext_loc(1)
                f_ext_y = f_ext_loc(2)
                f_ext_z = f_ext_loc(3)
              end where

            else if(compensate_f_ext) then ! force applied to a central particle only
                pd = input_int("dominika_particle_diameter",1)
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
                do i=n1/2+1-pdr,n1/2+1+pdr
                  do j=n2/2+1-pdr,n2/2+1+pdr
                    do k=n3/2+1-pdr,n3/2+1+pdr
                      if( norm2(real([ i-(n1/2+1), j-(n2/2+1), k-(n3/2+1) ],dp)) > real(pd,dp)/2._dp ) cycle
                      if (nature(i,j,k)/=fluid) err=.true.
                      f_ext_x(i,j,k) = f_ext_loc(1)
                      f_ext_y(i,j,k) = f_ext_loc(2)
                      f_ext_z(i,j,k) = f_ext_loc(3)
                      l=l+1
                      WRITE(47,*)i,j,k ! use ListPointPlot3D[data,BoxRatios->{1,1,1}] in Mathematica to read this file
                    end do
                  end do
                end do
                close(47)
                if(err.eqv..true.) then
                  print*,"ERROR: l306 of equilibration_new.f90. Dominika's particle at a solid node"
                  stop
                end if

                where(f_ext_x==f_ext_loc(1) .and. f_ext_y==f_ext_loc(2) .and. f_ext_z==f_ext_loc(3) )
                  f_ext_x = f_ext_x/l
                  f_ext_y = f_ext_y/l
                  f_ext_z = f_ext_z/l
                else where
                  f_ext_x = -f_ext_loc(1)/(fluid_nodes-l)
                  f_ext_y = -f_ext_loc(2)/(fluid_nodes-l)
                  f_ext_z = -f_ext_loc(3)/(fluid_nodes-l)
                end where

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
  DO j=1,ly
      DO k=1,lz
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
