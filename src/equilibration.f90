module module_equilibration
    implicit none
    private
    public :: equilibration
contains
    subroutine equilibration(n)

    USE precision_kinds, only: dp
    USE system, only: fluid, node, lbm, pbc
    use module_collision, only: collide
    use module_input, only: getinput
    USE mod_time, only: tick, tock

    implicit none
    real(dp), intent(inout), contiguous :: n(:,:,:,:) ! xyz;l
    integer :: t,i,j,k,l,ip,jp,kp,nx,ny,nz, lmin, lmax, pdr, pd, ios, px, py, pz
    integer :: nfluid, print_frequency, tfext, print_files_frequency, GL
    integer(kind(fluid)), allocatable, dimension(:,:,:) :: nature
    real(dp) :: n_loc, f_ext_loc(3), l2err, target_error
    real(dp), allocatable, dimension(:,:,:) :: density, jx, jy, jz, n_old, jx_old, jy_old, jz_old, f_ext_x, f_ext_y, f_ext_z
    integer, allocatable, dimension(:) :: cx, cy, cz
    logical :: convergence_reached, compensate_f_ext, convergence_reached_without_fext, convergence_reached_with_fext, err
    REAL(dp), PARAMETER :: eps=EPSILON(1._dp)
    LOGICAL :: write_total_mass_flux
    integer, allocatable :: il(:,:), jl(:,:), kl(:,:), l_inv(:)
    real(dp), parameter :: zerodp=0._dp
    real(dp) :: t0, t1, t2, t3, t4, t5

    ! Check system sizes.
    nx = getinput%int("lx", assert=">0")
    ny = getinput%int("ly", assert=">0")
    nz = getinput%int("lz", assert=">0")
    if( size(n,1)/=nx ) error stop "nx is inconsistant with n in equilibration"
    if( size(n,2)/=ny ) error stop "ny is inconsistant with n in equilibration"
    if( size(n,3)/=nz ) error stop "nz is inconsistant with n in equilibration"

    ! check LatticeBoltzmann model (we count velocity indices from 1, Fortran style!)
    lmin = lbm%lmin
    lmax = lbm%lmax
    if( lbound(n,4)/=1 .or. lmin/=1) error stop "lmin is strange and inconsistant in equilibration"
    if( ubound(n,4)/=lmax ) error stop "lmax is strange and inconsistant in equilibration"

    ! Print info to terminal every that number of steps
    print_frequency = getinput%int('print_frequency', defaultvalue=max(int(50000/(nx*ny*nz)),1), assert=">0" ) ! this number is my own optimal. To be generalized on strong criteria some day.

    ! Print velocity profiles to terminal every that number of steps
    print_files_frequency = getinput%int("print_files_frequency", HUGE(1))

    ! number of fluid nodes
    nfluid = count( node%nature==fluid )

    ! density
    allocate( density(nx,ny,nz), source=sum(n,4), stat=ios); if (ios /= 0) stop "density: Allocation request denied"

    ! l_inv gives you the inverse of the velocity v_l
    allocate( l_inv(lmin:lmax), source=lbm%vel(:)%inv, stat=ios); if(ios/=0) error stop "pb alloc l_inv in equilibration.f90"

    ! jx=ni*ci. j_old are used for checking convergence only
    allocate( jx     (nx,ny,nz), source=0._dp)
    allocate( jx_old (nx,ny,nz), source=0._dp)
    allocate( jy     (nx,ny,nz), source=0._dp)
    allocate( jy_old (nx,ny,nz), source=0._dp)
    allocate( jz     (nx,ny,nz), source=0._dp)
    allocate( jz_old (nx,ny,nz), source=0._dp)


    l2err = 999.
    open(13,file="./output/l2err.dat")

    OPEN(66, FILE="output/mass-flux_profile_along_z.dat")
    OPEN(67, FILE="output/mass-flux_profile_along_y.dat")
    OPEN(68, FILE="output/mass-flux_profile_along_x.dat")
    WRITE(66,*) "# z, <ρ.v_x>_{x,y}, <ρ.v_y>_{x,y}, <ρ.v_z>_{x,y}"
    WRITE(67,*) "# y, <ρ.v_x>_{x,z}, <ρ.v_y>_{x,z}, <ρ.v_z>_{x,z}"
    WRITE(68,*) "# x, <ρ.v_x>_{y,z}, <ρ.v_y>_{y,z}, <ρ.v_z>_{y,z}"
    OPEN(56, FILE="output/mean-density_profile_along_z.dat")
    OPEN(57, FILE="output/mean-density_profile_along_y.dat")
    OPEN(58, FILE="output/mean-density_profile_along_x.dat")
    write_total_mass_flux = getinput%log("write_total_mass_flux", .FALSE.)
    IF( write_total_mass_flux ) OPEN( 65, FILE="output/total_mass_flux.dat" )

    allocate( cx(lmax), source=lbm%cx )
    allocate( cy(lmax), source=lbm%cy )
    allocate( cz(lmax), source=lbm%cz )

    allocate( nature (nx,ny,nz), source=node%nature)
    allocate( f_ext_x(nx,ny,nz), source=zerodp)
    allocate( f_ext_y(nx,ny,nz), source=zerodp)
    allocate( f_ext_z(nx,ny,nz), source=zerodp)

    ! The force to be applied initialy: 0.
    f_ext_loc = zerodp

    !
    ! Tabulate the index of the node one finishes if one starts from a node and a velocity index l
    ! per direction
    !
    ! il is the index in x direction to which points a velocity cl from node i.
    allocate( il(lbm%lmin:lbm%lmax, 1:nx), stat=ios); if (ios /= 0) stop "il: Allocation request denied"
    allocate( jl(lbm%lmin:lbm%lmax, 1:ny), stat=ios); if (ios /= 0) stop "jl: Allocation request denied"
    allocate( kl(lbm%lmin:lbm%lmax, 1:nz), stat=ios); if (ios /= 0) stop "kl: Allocation request denied"
    do l= lmin, lmax
        il(l,:) = [( pbc(i+cx(l),1) ,i=1,nx )]
        jl(l,:) = [( pbc(j+cy(l),2) ,j=1,ny )]
        kl(l,:) = [( pbc(k+cz(l),3) ,k=1,nz )]
    end do

    convergence_reached_without_fext = .false.
    convergence_reached_with_fext = .false.

    !
    ! Compensate_f_ext is an option to have a local force one could apply on a given node (only the central node for now)
    ! compensated by a continuum background force, a little bit like a compensating electric field in a charged supercell.
    !
    ! Ade: By doing so, we can apply a force fx or fy in order to analyse and observe the velocity streamlines around
    ! the particle in the output file v_centralnode.dat.
    !
    compensate_f_ext = getinput%log("compensate_f_ext", defaultvalue=.false.)
    if(compensate_f_ext) call init_work_for_Adelchi()

    target_error = getinput%dp("target_error", 1.D-10)
    PRINT*
    PRINT*,'Lattice Boltzmann'
    PRINT*,'================='
    PRINT*,'       step     error   (target=',real(target_error),')'
    PRINT*,'       ---------------------------------------------'

    !
    ! TIME STEPS (in lb units)
    !
    do t=1,HUGE(t)
        call lots_of_prints()
        !
        ! Collision step
        !
CALL CPU_TIME(t0)
        call collide(n, jx, jy, jz, f_ext_x, f_ext_y, f_ext_z)
CALL CPU_TIME(t1)
        ! print velocity profile if you need/want it
        ! if( modulo(t, print_frequency) == 0) then
        !    call velocity_profiles(t) ! print velocity profiles
        ! end if

        !
        ! Bounce back (boundpm) to simplify propagation step
        !
        select case(GL)
        case(-1) ! bulk fluid, thus no bounce back
            ! do nothing
        case(1)! slit normal to z
            do l=lmin,lmax,2
                do k=1,nz
                    if( k/=1 .and. k/=nz) cycle
                    kp = kl(l,k)
                    !kp=pbc(k+cz(l),z)
                    do j=1,ny
                        jp = jl(l,j)
                        !jp=pbc(j+cy(l),y)
                        do i=1,nx
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
        case default ! generic case, we know nothing about the geometry a priori
            do l=lmin,lmax,2
                do k=1,nz
                    kp = kl(l,k)
                    !kp=pbc(k+cz(l),z)
                    do j=1,ny
                        jp = jl(l,j)
                        !jp=pbc(j+cy(l),y)
                        do i=1,nx
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
        end select

        !
        ! propagation
        !
        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(n,nz,ny,nx,lmin,lmax,il,jl,kl) &
        !$OMP PRIVATE(l,k,j,i,ip,jp,kp,n_old)
        do l=lmin,lmax
            n_old = n(:,:,:,l)
            do k=1,nz
                kp = kl(l,k)
                do j=1,ny
                    jp = jl(l,j)
                    do i=1,nx
                        ip = il(l,i)
                        n(ip,jp,kp,l) = n_old(i,j,k)
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO
        !
        ! The populations, that are probabilities, must never be negative
        !
        if( any(n<0) ) error stop "In equilibration_new, the population n(x,y,z,vel) < 0"

        IF( write_total_mass_flux ) WRITE(65,*) t, SUM(jx), SUM(jy), SUM(jz)

        !
        ! backup mass flux to test convergence at the end of the timestep
        !
        jx_old = jx
        jy_old = jy
        jz_old = jz

        ! update momentum densities after the propagation
        ! this is completely local in space and my be parallelized very well
        jx=f_ext_x/2.
        jy=f_ext_y/2.
        jz=f_ext_z/2.
        ! !$OMP PARALLEL DO DEFAULT(NONE)&
        ! !$OMP PRIVATE(l)&
        ! !$OMP SHARED(lmin,lmax,n,cx,cy,cz)&
        ! !$OMP REDUCTION(+:jx)&
        ! !$OMP REDUCTION(+:jy)&
        ! !$OMP REDUCTION(+:jz)
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
            do i=1,nx
                do k=1,nz
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
            do i=1,nx
                do k=1,nz
                    WRITE(90,*) i, k, f_ext_x(i,py,k), f_ext_z(i,py,k)
                    WRITE(91,*) i, k, jx(i,py,k), jz(i,py,k)
                end do
            end do
            close(90)
            close(91)
        end if

        !
        ! Check convergence
        ! Note to myself: we store these large arrays jx_old etc just for this!
        !
        l2err = maxval(   [  maxval(abs(jx-jx_old)), maxval(abs(jy-jy_old)), maxval(abs(jz-jz_old)) ])
        write(13,*) t, l2err
        if( l2err <= target_error .and. t>2 ) then
            convergence_reached = .true.
        else
            convergence_reached = .false.
        end if

        !
        ! Applying the forces or not?
        !
        if(convergence_reached) then
          if( .not.convergence_reached_without_fext ) then
            convergence_reached_without_fext = .true.
          else if( convergence_reached_without_fext ) then
            convergence_reached_with_fext = .true.
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
                f_ext_x = zerodp
                f_ext_y = zerodp
                f_ext_z = zerodp

                l=0 ! ADE: l counts the number of node within the particle
                err=.false.
                open(47,file="output/dominika_particle_shape.xyz")
                open(14,file="output/NodesInParticle.dat")
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
                  f_ext_x = -f_ext_loc(1)/(nfluid) +f_ext_x/l
                  f_ext_y = -f_ext_loc(2)/(nfluid) +f_ext_y/l
                  f_ext_z = -f_ext_loc(3)/(nfluid) +f_ext_z/l
                 else where
                  f_ext_x = -f_ext_loc(1)/(nfluid)
                  f_ext_y = -f_ext_loc(2)/(nfluid)
                  f_ext_z = -f_ext_loc(3)/(nfluid)
                 end where
                   if( any(abs([sum(f_ext_x)/nfluid,sum(f_ext_y)/nfluid,sum(f_ext_z)/nfluid])> eps ) ) then
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
  density=sum(n,4)
  WRITE(66,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(67,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(68,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(56,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(57,*)"# Steady state with convergence criteria", REAL(target_error)
  WRITE(58,*)"# Steady state with convergence criteria", REAL(target_error)
  DO k=1,nz
      WRITE(66,*) k, SUM(jx(:,:,k)), SUM(jy(:,:,k)), SUM(jz(:,:,k))
      WRITE(56,*) k, SUM(density(:,:,k))/ MAX( COUNT(density(:,:,k)>eps) ,1)
  END DO
  DO k=1,ny
      WRITE(67,*) k, SUM(jx(:,k,:)), SUM(jy(:,k,:)), SUM(jz(:,k,:))
      WRITE(57,*) k, SUM(density(:,k,:))/ MAX( COUNT(density(:,k,:)>eps) ,1)
  END DO
  DO k=1,nx
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
  DO j=1,ny
      DO k=1,nz
          WRITE(69,*) j, k, jy(1,j,k), jz(1,j,k)
      END DO
  END DO
  CLOSE(69)


  if( compensate_f_ext ) then
    print*,"       The particle is located at ", px,py,pz
    open(90,file="./output/f_ext-field.dat")
    open(91,file="./output/vel-field_central.dat")
    do i=1,nx
      do k=1,nz
        WRITE(90,*) i, k, f_ext_x(i,py,k), f_ext_z(i,py,k)
        WRITE(91,*) i, k,      jx(i,py,k),      jz(i,py,k)
      end do
    end do
    close(90)
    close(91)
  end if

  ! put back arrays into original types
  node%solventdensity = density
  node%solventflux(1) = jx
  node%solventflux(2) = jy
  node%solventflux(3) = jz
  node%nature = nature

contains

subroutine init_work_for_Adelchi()
  integer :: pCoord(3)
  ! We log the velocity of the particle node
  open(79,file="./output/v_centralnode.dat")
  pCoord = getinput%int3("particle_coordinates", defaultvalue=[nx/2+1,ny/2+1,nz/2+1], assert=">0" )
  px = pCoord(1)
  py = pCoord(2)
  pz = pCoord(3)
  if(px<1 .or. py<1 .or. pz<1 .or. px>nx .or. py>ny .or. pz>nz) error stop "Problem in central particle coordinates"
  pd = getinput%int("dominika_particle_diameter",1)
  print*,"Dominika's particle has a diameter of ", pd,"lb units"
  if( modulo(pd,2)==0 ) then
    print*,"ERROR: l. 285 particle diameter must be odd"
    print*,"-----  It is now",pd
    stop
  end if
  if(modulo(nx,2)==0 .or. modulo(ny,2)==0 .or. modulo(nz,2)==0) then
    print*,"ERROR: l.158 of equilibration_new.f90"
    print*,"=====  when compensate_f_ext, there should be odd number of nodes in all directions"
    print*,"nx, ny, nz =",nx,ny,nz
    stop
  end if
  pdr = pd/2 ! nodes of the particle on the right (or left) of the particle center. If particle is diameter 3, we have 1 node on the left and 1 on the right, so pd=3, pdr=3/2=1
end subroutine init_work_for_Adelchi


subroutine lots_of_prints
if( modulo(t, print_frequency) == 0) PRINT*, t, l2err
if( MODULO(t, print_files_frequency) == 0 .OR. t==1) THEN
    WRITE(66,*)"# timestep",t
    WRITE(67,*)"# timestep",t
    WRITE(68,*)"# timestep",t
    WRITE(56,*)"# timestep",t
    WRITE(57,*)"# timestep",t
    WRITE(58,*)"# timestep",t
    DO k=1,nz
        WRITE(56,*) k, SUM(density(:,:,k))/ MAX( COUNT(density(:,:,k)>eps)  ,1)
        WRITE(66,*) k, SUM(jx(:,:,k)), SUM(jy(:,:,k)), SUM(jz(:,:,k))
    END DO
    DO k=1,ny
        WRITE(57,*) k, SUM(density(:,k,:))/ MAX( COUNT(density(:,k,:)>eps)  ,1)
        WRITE(67,*) k, SUM(jx(:,k,:)), SUM(jy(:,k,:)), SUM(jz(:,k,:))
    END DO
    DO k=1,nx
        WRITE(58,*) k, SUM(density(k,:,:))/ MAX( COUNT(density(k,:,:)>eps)  ,1)
        WRITE(68,*) k, SUM(jx(k,:,:)), SUM(jy(k,:,:)), SUM(jz(k,:,:))
    END DO
    WRITE(66,*)
    WRITE(67,*)
    WRITE(68,*)
END IF

if( compensate_f_ext .and. convergence_reached_without_fext) then
    write(79,*)t-tfext, jx(px,py,pz), jy(px,py,pz), jz(px,py,pz)
end if
end subroutine lots_of_prints


end subroutine equilibration
end module module_equilibration
