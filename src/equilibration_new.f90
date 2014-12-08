subroutine equilibration_new
  use precision_kinds, only: i2b, dp, sp
  use system, only: fluid, supercell, node, lbm, n, pbc
  use populations, only: update_populations
  use input, only: input_dp3, input_dp, input_int
  use constants, only: x, y, z, zerodp, epsdp
  implicit none
  integer :: t,i,j,k,l,ip,jp,kp,n1,n2,n3,lmax,tmax,l_inv
  integer :: fluid_nodes, print_frequency
  integer(kind(fluid)), allocatable, dimension(:,:,:) :: nature
  real(dp) :: sigma, n_loc, f_ext(3), f_ext_loc(3), l2err, target_error
  real(dp), allocatable, dimension(:,:,:) :: density, jx, jy, jz, old_n, jx_old, jy_old, jz_old
  real(dp), allocatable, dimension(:) :: a0, a1
  integer, allocatable, dimension(:) :: cx, cy, cz
  logical :: convergence_reached

  sigma = input_dp('sigma', zerodp) ! net charge of the solid phase. Kind of an external potential.
  if( abs(sigma) > epsilon(1._dp) ) then
    print*,"ERROR: laboetie is valid only for uncharged systems."
    print*,"===== Dont tell Benjamin you'd like to see such feature in Laboetie :)"
    print*,"Hi Benjamin ! :))"
    stop
  end if

  print_frequency = input_int('print_frequency',10000)
  fluid_nodes = count( node%nature==fluid )
  f_ext = zerodp ! external potential (pressure gradient init to 0)
  tmax = input_int("tmax") ! maximum number of iterations
  target_error = input_dp("target_error")
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
  allocate( nature (n1,n2,n3), source=node%nature)
  lmax = lbm%lmax
  allocate( cx(lmax), source=lbm%vel(:)%coo(1))
  allocate( cy(lmax), source=lbm%vel(:)%coo(2))
  allocate( cz(lmax), source=lbm%vel(:)%coo(3))
  allocate( a0(lmax), source=lbm%vel(:)%a0)
  allocate( a1(lmax), source=lbm%vel(:)%a1)

  print*,'       step     <j>_x            <j>_y            <j>_z             err.            target.err.'
  print*,'       ----------------------------------------------------------------------------------------'

  do t=1,huge(t) !tmax

    if( modulo(t, print_frequency) == 0) then
      print*,t,&
      real(sum(jx/density, mask=node%nature==fluid)/fluid_nodes,sp), &
      real(sum(jy/density, mask=node%nature==fluid)/fluid_nodes,sp), &
      real(sum(jz/density, mask=node%nature==fluid)/fluid_nodes,sp), &
      real(l2err,sp), real(target_error,sp)
    end if

    ! backup moment density (velocities) to test convergence at the end of the timestep
    jx_old = jx
    jy_old = jy
    jz_old = jz

    ! collision
    do concurrent(i=1:n1, j=1:n2, k=1:n3, l=1:lmax)
      if( node(i,j,k)%nature == fluid ) then
        f_ext_loc = f_ext
      else
        f_ext_loc = zerodp
      end if
      n(i,j,k,l) = a0(l)*density(i,j,k) +a1(l)*(&
      cx(l)*(jx(i,j,k)+f_ext_loc(x)) + cy(l)*(jy(i,j,k)+f_ext_loc(y)) + cz(l)*(jz(i,j,k)+f_ext_loc(z)))
      ! n(i,j,k,l) = a0(l)*density(i,j,k) +a1(l)* dot_product( &
        ! [cx(l), cy(l), cz(l)],&
        ! [jx(i,j,k)+f_ext_loc(x) , jy(i,j,k)+f_ext_loc(y) , jz(i,j,k)+f_ext_loc(z)] )
    end do

    ! print velocity profile if you need/want it
    ! if( modulo(t, print_frequency) == 0) then
    !    call velocity_profiles(t) ! print velocity profiles
    ! end if

    ! bounce back (boundpm) to simplify propagation step
    do concurrent( i=1:n1, j=1:n2, k=1:n3, l=1:lmax:2)
      ip=pbc(i+cx(l),x)
      jp=pbc(j+cy(l),y)
      kp=pbc(k+cz(l),z)
      if( nature(i,j,k) /= nature(ip,jp,kp) ) then
        l_inv = lbm%vel(l)%inv
        n_loc = n(i,j,k,l)
        n(i,j,k,l) = n(ip,jp,kp,l_inv)
        n(ip,jp,kp,l_inv) = n_loc
      end if
    end do

    ! propagation step
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

    ! check new populations
    if(any(n<0)) then
      print*,"ERROR: n(i,j,k,l) is negative somewhere. Check l.115 of equilibration_new.f90"
      print*,"====="
      stop
    end if

    ! update densities after the propagation and check it
    ! this is also completely local in space
    density = sum(n,4)

    ! print*,sum(density)
    ! if( abs(sum(density)/real(n1*n2*n3,kind(density)) -1._dp) > epsilon(1._dp) ) then
    !   stop "otto"
    ! end if

    ! update momentum densities after the propagation
    ! this is completely local in space and my be parallelized very well
    do concurrent (i=1:n1, j=1:n2, k=1:n3)
      jx(i,j,k) = (jx(i,j,k) + sum(n(i,j,k,:)*cx(:)))/2._dp
      jy(i,j,k) = (jy(i,j,k) + sum(n(i,j,k,:)*cy(:)))/2._dp
      jz(i,j,k) = (jz(i,j,k) + sum(n(i,j,k,:)*cz(:)))/2._dp
    end do

    ! check convergence
    l2err = norm2(jx-jx_old+jy-jy_old+jz-jz_old)
    if( l2err <= target_error ) then
      convergence_reached = .true.
    else
      convergence_reached = .false.
    end if

    ! chose to apply external contraints (f_ext) or not
    if( convergence_reached ) then
      if( norm2(f_ext-input_dp3("f_ext")) <= epsdp .and. t>2) then ! if equilibration with constraints is already converged
        exit
      else ! if equilibration is without external forces, then now add external forces
        print*,"       Applying constraints"
        f_ext = input_dp3("f_ext")
      end if
    end if

  end do

  print*,"       Convergence reached at time step",t-1

end subroutine equilibration_new
