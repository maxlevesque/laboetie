subroutine equilibration_new
  use precision_kinds, only: i2b, dp, sp
  use system, only: fluid, supercell, node, lbm, n, pbc
  use populations, only: update_populations
  use input, only: input_dp3, input_dp, input_int, input_log
  use constants, only: x, y, z, zerodp, epsdp
  implicit none
  integer :: t,i,j,k,l,ip,jp,kp,n1,n2,n3,lmax,tmax,l_inv
  integer :: fluid_nodes, print_frequency, supercellgeometrylabel, tfext
  integer(kind(fluid)), allocatable, dimension(:,:,:) :: nature
  real(dp) :: sigma, n_loc, f_ext_loc(3), l2err, target_error, minimumvalueofJ, maximumvalueofJ
  real(dp), allocatable, dimension(:,:,:) :: density, jx, jy, jz, old_n, jx_old, jy_old, jz_old, f_ext_x, f_ext_y, f_ext_z
  real(dp), allocatable, dimension(:) :: a0, a1
  integer, allocatable, dimension(:) :: cx, cy, cz
  logical :: convergence_reached, compensate_f_ext, convergence_reached_without_fext, convergence_reached_with_fext

  sigma = input_dp('sigma', zerodp) ! net charge of the solid phase. Kind of an external potential.
  if( abs(sigma) > epsilon(1._dp) ) then
    print*,"ERROR: laboetie is valid only for uncharged systems."
    print*,"===== Dont tell Benjamin you'd like to see such feature in Laboetie :)"
    print*,"Hi Benjamin ! :))"
    stop
  end if

  supercellgeometrylabel = supercell%geometry%label ! -1 for solid free cell

  print_frequency = input_int('print_frequency',10000)
  fluid_nodes = count( node%nature==fluid )
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
  allocate( f_ext_x(n1,n2,n3), source=zerodp)
  allocate( f_ext_y(n1,n2,n3), source=zerodp)
  allocate( f_ext_z(n1,n2,n3), source=zerodp)
  f_ext_loc = zerodp ! this is important and spagetty like... please read carefully before modifying this line
  lmax = lbm%lmax
  allocate( cx(lmax), source=lbm%vel(:)%coo(1))
  allocate( cy(lmax), source=lbm%vel(:)%coo(2))
  allocate( cz(lmax), source=lbm%vel(:)%coo(3))
  allocate( a0(lmax), source=lbm%vel(:)%a0)
  allocate( a1(lmax), source=lbm%vel(:)%a1)

  convergence_reached_without_fext = .false.
  convergence_reached_with_fext = .false.
  compensate_f_ext = input_log("compensate_f_ext")
  if(compensate_f_ext) open(79,file="./output/v_centralnode.dat")

  print*,'       step    maxval(j)          L2.err.          target.err.'
  print*,'       ----------------------------------------------------------------------------------------'

  ! TIME STEPS
  do t=1,huge(t)

    if( modulo(t, print_frequency) == 0) then
      minimumvalueofJ = min(  minval(abs(jx)/density,density>epsdp),&
                              minval(abs(jy)/density,density>epsdp),&
                              minval(abs(jz)/density,density>epsdp)  )
      if(minimumvalueofJ<=epsdp) minimumvalueofJ = zerodp
      maximumvalueofJ = max(  maxval(abs(jx)/density,density>epsdp),&
                              maxval(abs(jy)/density,density>epsdp),&
                              maxval(abs(jz)/density,density>epsdp)  )
      print*, t, real([minimumvalueofJ,maximumvalueofJ,l2err,target_error],sp)
    end if

    ! VACF of central node
    if( compensate_f_ext .and. convergence_reached_without_fext) then
      write(79,*)t-tfext, jz(n1/2+1,n2/2+1,n3/2+1)
    end if

    ! backup moment density (velocities) to test convergence at the end of the timestep
    jx_old = jx
    jy_old = jy
    jz_old = jz

    ! collision
    do concurrent(l=1:lmax)
      n(:,:,:,l) = a0(l)*density +a1(l)*(cx(l)*(jx+f_ext_x)+cy(l)*(jy+f_ext_y)+cz(l)*(jz+f_ext_z))
    end do
    ! do concurrent(i=1:n1, j=1:n2, k=1:n3, l=1:lmax)
    !   n(i,j,k,l) = a0(l)*density(i,j,k) +a1(l)*(&
    !     cx(l)*(jx(i,j,k)+f_ext_x(i,j,k)) + cy(l)*(jy(i,j,k)+f_ext_y(i,j,k)) + cz(l)*(jz(i,j,k)+f_ext_z(i,j,k)))
    ! end do

    ! print velocity profile if you need/want it
    ! if( modulo(t, print_frequency) == 0) then
    !    call velocity_profiles(t) ! print velocity profiles
    ! end if

    ! BOUNCE BACK (boundpm) to simplify propagation step
    if(supercellgeometrylabel/=-1) then ! if supercell has fluid nodes only, bounce back is useless
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
    end if

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
    do concurrent (l=1:lmax)
      jx = jx +n(:,:,:,l)*cx(l)
      jy = jy +n(:,:,:,l)*cy(l)
      jz = jz +n(:,:,:,l)*cz(l)
    end do
    jx=jx/2
    jy=jy/2
    jz=jz/2
    ! do concurrent (i=1:n1, j=1:n2, k=1:n3)
    !   jx(i,j,k) = (jx(i,j,k) + sum(n(i,j,k,:)*cx(:)))/2._dp
    !   jy(i,j,k) = (jy(i,j,k) + sum(n(i,j,k,:)*cy(:)))/2._dp
    !   jz(i,j,k) = (jz(i,j,k) + sum(n(i,j,k,:)*cz(:)))/2._dp
    ! end do



    if( compensate_f_ext .and. convergence_reached_without_fext .and. t==tfext) then
      open(90,file="./output/f_ext-field_t0.dat")
      open(91,file="./output/vel-field_central_t0.dat")
      do i=1,n1
        do k=1,n3
          write(90,*) i, k, f_ext_x(i,n2/2+1,k), f_ext_z(i,n2/2+1,k)
          write(91,*) i, k,      jx(i,n2/2+1,k)     ,      jz(i,n2/2+1,k)
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
          write(90,*) i, k, f_ext_x(i,n2/2+1,k), f_ext_z(i,n2/2+1,k)
          write(91,*) i, k, jx(i,n2/2+1,k), jz(i,n2/2+1,k)
        end do
      end do
      close(90)
      close(91)
    end if





    ! check convergence
    l2err = norm2(jx-jx_old+jy-jy_old+jz-jz_old)
    if( l2err <= target_error .and. t>2) then
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

    ! chose to apply external contraints (f_ext) or not
    if( convergence_reached ) then

      ! if you are already converged without then with f_ext then quit time loop. Stationary state is found.
      if( convergence_reached_without_fext .and. convergence_reached_with_fext .and. t>2) then
        exit ! loop over time steps

      ! if you have already converged without fext, but not yet with fext, then enable fext
      else if(convergence_reached_without_fext .and. .not.convergence_reached_with_fext) then
        tfext=t+1
        print*,"       Applying constraints at time step",tfext
        f_ext_loc = input_dp3("f_ext")
        print*,"       Pressure gradient (f_ext in lb.in) is",f_ext_loc
        if(.not.compensate_f_ext) then ! the force is exerced everywhere with same intensity
          print*,"       It is applied homogeneously everywhere in the fluid"
          f_ext_x = f_ext_loc(1)
          f_ext_y = f_ext_loc(2)
          f_ext_z = f_ext_loc(3)
        else if(compensate_f_ext) then ! force applied to central node only
          print*,"       It is applied to the central node only. Other nodes see a homogeneous compensating potential"
          if(modulo(n1,2)==0 .or. modulo(n2,2)==0 .or. modulo(n3,2)==0) then
            print*,"ERROR: l.158 of equilibration_new.f90"
            print*,"=====  when compensate_f_ext, there should be odd number of nodes in all directions"
            print*,"n1, n2, n3 =",n1,n2,n3
            stop
          end if
          where(nature==fluid)
            f_ext_x = -f_ext_loc(1)/(fluid_nodes-1)
            f_ext_y = -f_ext_loc(2)/(fluid_nodes-1)
            f_ext_z = -f_ext_loc(3)/(fluid_nodes-1)
          else where
            f_ext_x = zerodp
            f_ext_y = zerodp
            f_ext_z = zerodp
          end where
          if(nature(n1/2+1,n2/2+1,n3/2+1)/=fluid) then
            print*,"ERROR: l.168 of equilibration_new.f90"
            print*,"=====  The central node must be fluid. Otherwise is madness! :)"
            print*,"This is sparta!"
            stop
          end if
          f_ext_x(n1/2+1,n2/2+1,n3/2+1) = f_ext_loc(1)
          f_ext_y(n1/2+1,n2/2+1,n3/2+1) = f_ext_loc(2)
          f_ext_z(n1/2+1,n2/2+1,n3/2+1) = f_ext_loc(3)
          if( any( abs([sum(f_ext_x),sum(f_ext_y),sum(f_ext_z)]) > epsdp ) ) then
            print*,"ERROR: l.215 of equilibration_new.f90"
            print*,"=====  The compensation is not well-implemented."
            stop
          end if
        end if
      end if
    end if

  end do

  close(79)
  print*,"       Convergence reached at time step",t-1

  if( compensate_f_ext ) then
    open(90,file="./output/f_ext-field.dat")
    open(91,file="./output/vel-field_central.dat")
    do i=1,n1
      do k=1,n3
        write(90,*) i, k, f_ext_x(i,n2/2+1,k), f_ext_z(i,n2/2+1,k)
        write(91,*) i, k,      jx(i,n2/2+1,k),      jz(i,n2/2+1,k)
      end do
    end do
    close(90)
    close(91)
  end if

end subroutine equilibration_new
