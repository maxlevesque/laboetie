MODULE MOMENT_PROPAGATION

  USE precision_kinds
  USE constants, only: x, y, z, zerodp
  USE system, only: tmax, tmin, tmom, pbc, supercell, delta_disp, delta_vacf, delta_fa, fluid, solid
  USE mod_lbmodel, only: lbm

  IMPLICIT NONE

  !PRIVATE
  !PUBLIC :: propagate, deallocate_propagated_quantity, Adsorption_equilibration 
  !PUBLIC :: init_adsorption_equilibration, init_propagation

  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: Propagated_Quantity, Propagated_Quantity_Adsorbed, Adsorbed_density, Free_density 
  INTEGER(i2b), PARAMETER :: now=0, next=1, past=-1, tini=past
  REAL(dp), DIMENSION(x:z, past:next) :: vacf
  REAL(dp) :: lambda, lambda_s, alpha
  TYPE type_tracer
    REAL(dp) :: ka, kd, K, z, Db, Ds, pmax !K=ka/kd, z=tracer charge
  END TYPE
  TYPE(type_tracer) :: tracer
  INTEGER(i2b), PRIVATE :: lx,ly,lz
  INTEGER(i2b) :: considerpmax, nf
  logical :: considerAdsorption
  REAL(dp), DIMENSION(3) :: disp_now, disp_past, deltad, diff_past, diff_now
  DOUBLE PRECISION :: fa_past, fa_now, deltaf, cs
  integer(i2b) :: i, j, k, l, ip, jp, kp, nature_loc, ll, lu, l_inv_loc
  integer(i2b), allocatable :: ip_all(:), jp_all(:), kp_all(:)
  integer(i2b), save, allocatable :: c(:,:), l_inv(:)
  integer(kind(fluid)), save, allocatable :: nature(:,:,:)
  real(dp), save, allocatable :: density(:,:,:), a0(:), a1(:)
  logical, save, allocatable :: interfacial(:,:,:)
  double precision, allocatable, dimension(:,:,:) :: surf_weight
  

  CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ADSORPTION INITIALIZATION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SUBROUTINE init_adsorption_equilibration

      use system, ONLY: fluid, solid, node
      use input, ONLY: input_dp, input_log, input_int
      implicit none
      integer(i2b) :: nf
      double precision :: spwx, spwy, spwz

      cs = 1.0d0/sqrt(3.0d0)
      ll = lbm%lmin
      lu = lbm%lmax
      lx = supercell%geometry%dimensions%indiceMax(x)
      ly = supercell%geometry%dimensions%indiceMax(y)
      lz = supercell%geometry%dimensions%indiceMax(z)
      if (.not. allocated(c)) then
        allocate(c(x:z,ll:lu))
        c(x,:) = lbm%vel(:)%coo(x)
        c(y,:) = lbm%vel(:)%coo(y)
        c(z,:) = lbm%vel(:)%coo(z)
      end if
      if (.not. allocated(a0)) allocate(a0(ll:lu), source=lbm%vel(:)%a0)
      if (.not. allocated(a1)) allocate(a1(ll:lu), source=lbm%vel(:)%a1)
      if (.not. allocated(l_inv)) allocate(l_inv(ll:lu), source=lbm%vel(:)%inv)
      if (.not. allocated(density)) allocate (density(lx,ly,lz), source=node(:,:,:)%solventDensity)
      if (.not. allocated(nature)) allocate (nature(lx,ly,lz), source=node(:,:,:)%nature)
      if (.not. allocated(interfacial)) allocate (interfacial(lx,ly,lz), source=node(:,:,:)%isInterfacial)
      if (.not. allocated(surf_weight)) allocate(surf_weight(lx,ly,lz), source=zerodp)
      if (.not. allocated(ip_all)) allocate(ip_all(lu))
      if (.not. allocated(jp_all)) allocate(jp_all(lu))
      if (.not. allocated(kp_all)) allocate(kp_all(lu))

      tracer%ka = input_dp('tracer_ka') ! adsorption
      tracer%kd = input_dp('tracer_kd') ! desorption
      if (tracer%ka < -epsilon(1._dp)) stop 'I detected tracer%ka to be <0 in module moment_propagation. STOP.'
      if (tracer%kd < -epsilon(1._dp)) stop 'I detected tracer%kd to be <0 in module moment_propagation. STOP.'

      if ( abs(tracer%kd)<=epsilon(1._dp) ) then ! if kd=0
        tracer%K = 0.0_dp ! don't do the division by 0
      else
        tracer%K = tracer%ka / tracer%kd
      end if
      
      if ( abs(tracer%K)>epsilon(1.0_dp) ) then
        considerAdsorption = .true.
      else
        considerAdsorption = .false.
      end if

      tracer%z = input_dp('tracer_z') ! tracer's charge
      tracer%Db = input_dp('tracer_Db') ! bulk diffusion coefficient of tracer, i.e. the molecular diffusion coefficient
      tracer%Ds = input_dp('tracer_Ds') ! surface diffusion coefficient of tracer
      IF (tracer%Db <= 0.0_dp ) STOP 'tracer_Db as readen in input is invalid'
      IF (tracer%Ds /= 0.0_dp ) STOP "I've found a non-zero Ds (surface diffusion coefficient) in input file. Not implemented yet"

      tracer%pmax = input_dp('pmax')
      considerpmax = input_int('considerpmax')
      if (considerpmax == 1) then
        alpha = 1.0d0/tracer%pmax
      else if (considerpmax == 0) then
       alpha = 0.0d0
      else
        STOP 'Pmax consideration is invalid!!'
      end if

      lambda = calc_lambda(tracer%Db)      ! bulk diffusion
      lambda_s = calc_lambda(tracer%Ds)  ! surface diffusion. is 0 for Ds=0
   
      call test_and_allocate_what_is_needed_for_moment_propagation
      
      nf = COUNT( node%nature==fluid )

       where (node%nature==fluid)
		 Free_density(1, :, :, :, now) = 1.0d0/(real(nf)*3.0d0)
		 Free_density(2, :, :, :, now) = 1.0d0/(real(nf)*3.0d0)
		 Free_density(3, :, :, :, now) = 1.0d0/(real(nf)*3.0d0)
	  end where
      print*, '1 - Fa = ', sum(Free_density(:,:,:,:,now))

      !!$OMP PARALLEL DO PRIVATE(i,j,k,l,ip,jp,kp,spwx,spwy,spwz)&
      !!$OMP PRIVATE(ip_all,jp_all,kp_all) &
      !!$OMP SHARED(nature,a0,c,ll,lu,lx,ly,lz, cs, surf_weight) &
      !!$OMP DEFAULT(NONE)
      do k=1,lz ! we parallelize over k. If system is 30x30x1 parallelization is useless!
        kp_all(:) = [( pbc( k+c(z,l) ,z) , l=ll, lu)]
        do j=1,ly
          jp_all(:) = [( pbc( j+c(y,l) ,y) , l=ll, lu)]
          do i=1,lx      
            ip_all(:) = [( pbc( i+c(x,l) ,x) , l=ll, lu)]  
            spwx = 0.0d0
            spwy = 0.0d0
            spwz = 0.0d0
            do l = ll+1, lu 
              ip = ip_all(l)
              jp = jp_all(l)
              kp = kp_all(l)
		      if (nature(i,j,k) == solid .and. nature(ip,jp,kp) == fluid) then 	    
                spwx = spwx + 2.0d0*a1(l)*c(x,l)
			    spwy = spwy + 2.0d0*a1(l)*c(y,l)
           	    spwz = spwz + 2.0d0*a1(l)*c(z,l)
              end if
            end do
            surf_weight(i,j,k) = sqrt(spwx**2 + spwy**2 + spwz**2)
          end do
        end do
      end do
      !!$OMP END PARALLEL DO
      print*, 'Surf_spe = ', sum(surf_weight)
      call plot_surf_weight

      PRINT*, ""
      print*, 'ADSORPTION EQUILIBRATION LOOP'
      print*, 'Time step      Free quantity         Adsorbed quantity       Total      Time(d:h:m:s)'

      OPEN(99, file='output/vacf.dat')
      OPEN(102, FILE='output/diff.dat')
      OPEN(103, FILE='output/disp.dat')
      OPEN(104, FILE='output/Fa.dat')
        WRITE(99,*) 0, vacf(x,past), vacf(y,past), vacf(z,past)
        WRITE(102,*) 0, 0.0, 0.0, 0.0
        WRITE(103,*) 0, 0.0, 0.0, 0.0
        WRITE(104,*) 0, 0.0, 0.0, 0.0
      CLOSE(99)
      CLOSE(102)
      CLOSE(103)
      CLOSE(104)
    END SUBROUTINE init_adsorption_equilibration

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!ADSORPTION EQUILIBRATION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Subroutine Adsorption_equilibration(it, is_converged)

      use system, only: fluid, solid, n, node, tic
      use input, only: input_char, input_log
      use OMP_lib
      implicit none
      integer(i2b), intent(in) :: it
      real(dp) :: fermi, fractionOfParticleRemaining, scattprop, scattprop_p, n_loc(lbm%lmin:lbm%lmax)
      real(dp) :: u_star(x:z), Free_density_loc(x:z), density_loc
      logical :: error, interfacial_loc
      logical, intent(out) :: is_converged
      integer(i4b) :: tac_loc, time_loc, count_rate, days, hours, minutes, seconds
      character(len=12) :: time_str, temp_str(4)

      error=.false.

      !$OMP PARALLEL DO PRIVATE(i,j,k,l,ip,jp,kp,scattprop,l_inv_loc,scattprop_p,n_loc,density_loc,nature_loc,u_star)&
      !$OMP PRIVATE(ip_all,jp_all,kp_all,fractionOfParticleRemaining, Free_density_loc, interfacial_loc, fermi) &
      !$OMP SHARED(nature,n,lambda, Free_density,l_inv,a0,c,ll,lu,lx,ly,lz,considerAdsorption,tracer) &
      !$OMP SHARED(density,Adsorbed_density,interfacial,error,alpha) &
      !$OMP REDUCTION(+:vacf) &
      !$OMP DEFAULT(NONE)
      
      

      do k=1,lz ! we parallelize over k. If system is 30x30x1 parallelization is useless!
        kp_all(:) = [( pbc( k+c(z,l) ,z) , l=ll, lu)]
        do j=1,ly
          jp_all(:) = [( pbc( j+c(y,l) ,y) , l=ll, lu)]
          do i=1,lx
            nature_loc = nature(i,j,k)
            if (nature_loc/=fluid) cycle
            ip_all(:) = [( pbc( i+c(x,l) ,x) , l=ll, lu)]
            u_star(:) = 0.0_dp ! the average velocity at r
            fractionOfParticleRemaining = 1.0_dp ! fraction of particles staying at r; decreases in the loop over neighbours
            n_loc(:) = n(i,j,k,:)
            density_loc = density(i,j,k)
            interfacial_loc = Interfacial(i,j,k)
            Free_density_loc(:) = Free_density(:,i,j,k,next)

            do l = ll+1, lu ! ll is velocity=0
              ip = ip_all(l)
              jp = jp_all(l)
              kp = kp_all(l)
              if ( nature(ip,jp,kp) /= fluid ) cycle
              fermi = 1.0_dp/(1.0_dp + calc_exp_dphi(i,j,k,ip,jp,kp))
              scattprop = calc_scattprop( n_loc(l), density_loc, a0(l), lambda, fermi) ! scattering probability at r
              fractionOfParticleRemaining = fractionOfParticleRemaining - scattprop ! what is scattered away is not found anymore at r
              u_star(:) = u_star(:) + scattprop*c(:,l)
              l_inv_loc = l_inv(l)
              scattprop_p = calc_scattprop( n(ip,jp,kp,l_inv_loc), density(ip,jp,kp), a0(l_inv_loc), lambda, 1.0_dp-fermi)
              Free_density_loc(:) = Free_density_loc(:) + Free_density(:,ip,jp,kp,now)*scattprop_p
            end do
            Free_density(:,i,j,k,next) = Free_density_loc(:)

            if (   (.not.interfacial_loc .and. considerAdsorption) &
              .or. (.not.considerAdsorption) )then
              Free_density(:,i,j,k,next) = &
                Free_density(:,i,j,k,next) + fractionOfParticleRemaining*Free_density(:,i,j,k,now)
            else if ( interfacial_loc .and. considerAdsorption ) then
              fractionOfParticleRemaining = fractionOfParticleRemaining&
                                            & - tracer%ka*(1 - sum(Adsorbed_density(:,i,j,k,now))*alpha) 
              Free_density(:,i,j,k,next) = Free_density(:,i,j,k,next) &
                + fractionOfParticleRemaining * Free_density(:,i,j,k,now) &
                + Adsorbed_density(:,i,j,k,now) * tracer%kd
              Adsorbed_density(:,i,j,k,next) = &
                Adsorbed_density(:,i,j,k,now) * (1.0_dp - tracer%kd) &
                + (1 - sum(Adsorbed_density(:,i,j,k,now))*alpha)*Free_density(:,i,j,k,now)*tracer%ka
            else
              stop "OMG YOU KILLED THE QUEEN!"
            end if

            if (abs(fractionOfParticleRemaining)<epsilon(1._dp)) error=.true.
          end do
        end do
      end do

!$OMP END PARALLEL DO
      
      if(error) stop 'somewhere restpart is negative' ! TODO one should find a better function for ads and des, as did Benjamin for pi

      Free_density(:,:,:,:,now) = Free_density(:,:,:,:,next)
      Free_density(:,:,:,:,next) = 0.0_dp
      if (considerAdsorption) then
        Adsorbed_density(:,:,:,:,now) = Adsorbed_density(:,:,:,:,next)
        Adsorbed_density(:,:,:,:,next) = 0.0_dp
      end if

      fa_now = sum(Adsorbed_density(:,:,:,:,now))
      deltaf = fa_now - fa_past
 
      if(modulo(it,10000)==0) then 
        CALL SYSTEM_CLOCK(tac_loc, count_rate)
        time_loc = (tac_loc-tic)/count_rate
        days = floor(time_loc/24.0/3600.0)
        time_loc = time_loc - days*24.0*3600
        hours = floor(time_loc/3600.0)
        time_loc = time_loc - hours*3600.0
        minutes = floor(time_loc/60.0)
        seconds = time_loc - minutes*60.0
        write(temp_str(1),'(I2.2)') days
        write(temp_str(2),'(I2.2)') hours
        write(temp_str(3),'(I2.2)') minutes
        write(temp_str(4),'(I2.2)') seconds
        time_str = trim(temp_str(1))//':'//trim(temp_str(2))//':'//trim(temp_str(3))//':'//trim(temp_str(4))
        !time_str = 'coucou'
        OPEN(104, FILE='output/Fa.dat', ACCESS='append')
            print*, it,&
                    real(sum(Free_density(:,:,:,:,now))),&
                    real(Fa_now),&
                    real(sum(Free_density(:,:,:,:,now)) + sum(Adsorbed_density(:,:,:,:,now))),&
                    time_str   
            WRITE(104,*) it,&
                    real(sum(Free_density(:,:,:,:,now))),&
                    real(Fa_now),&
                    real(sum(Free_density(:,:,:,:,now)) + sum(Adsorbed_density(:,:,:,:,now))),&
                    time_str
        CLOSE(104)
    end if

      fa_past = fa_now
    
    IF (it> 10001 .and. deltaf < delta_fa) THEN
        print*, 'CONVERGENCE IS REACHED: delta_fa <', delta_fa
        print*, ""
        print*, 'Fa = ', Fa_now
        is_converged = .true.
        call plot_densities_of_particles(it)
    ELSE 
        is_converged = .false.
        if (it == tmax) then
            call plot_densities_of_particles(it)
            print*, '!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!'
            print*, '!!!!!!!!!NO CONVERGENCE!!!!!!!!!!!'
            print*, ""
            print*, 'Fa = ', Fa_now
        end if    
    END IF

    end subroutine Adsorption_equilibration

SUBROUTINE init_propagation

      use system, ONLY: phi, fluid, solid, n, node
      use input, ONLY: input_dp, input_log, input_int
      implicit none
      real(dp) :: boltz_weight, scattprop, scattprop_p, fermi, exp_dphi, exp_min_dphi, rho
      real(dp) :: n_loc(lbm%lmin:lbm%lmax)
      integer(i2b) :: i, j, k, l, l_inv, ip, jp, kp
        
      diff_past = 0.0_dp
      diff_now = 0.0_dp
      disp_past = 0.0_dp
      disp_now = 0.0_dp
      deltad = 3000
      deltaf = 3000
      fa_now = 0.0_dp
      fa_past = 0.0_dp

      do concurrent (i=1:lx, j=1:ly, k=1:lz, node(i,j,k)%nature==fluid )
        boltz_weight = exp(-tracer%z*phi(i,j,k))*sum(Free_density(:,i,j,k,now))
        n_loc = n(i,j,k,:)
        rho = node(i,j,k)%solventDensity

        do concurrent (l=lbm%lmin+1:lbm%lmax)
          ip = pbc (i+lbm%vel(l)%coo(x) ,x)
          jp = pbc (j+lbm%vel(l)%coo(y) ,y)
          kp = pbc (k+lbm%vel(l)%coo(z) ,z)

          if (node(ip,jp,kp)%nature==solid) cycle
          exp_dphi = calc_exp_dphi( i, j, k, ip, jp, kp)
          exp_min_dphi = 1.0_dp/exp_dphi ! =1 if tracer%z=0
          fermi = 1.0_dp/(1.0_dp + exp_dphi) ! =0.5 if tracer%z=0
          scattprop = calc_scattprop( n_loc(l), rho, lbm%vel(l)%a0, lambda, fermi)

          vacf(:,tini) = vacf(:,tini) + boltz_weight*scattprop *lbm%vel(l)%coo(:)**2

          l_inv = lbm%vel(l)%inv ! comes to r
           scattprop_p = calc_scattprop( &
            n(ip,jp,kp,l_inv), node(ip,jp,kp)%solventDensity, lbm%vel(l_inv)%a0, lambda, 1.0_dp-fermi)

          Propagated_Quantity(:,i,j,k,tini+1) = Propagated_Quantity(:,i,j,k,tini+1) &
            + exp_min_dphi * scattprop_p * lbm%vel(l_inv)%coo(:) * boltz_weight

        end do
      end do
      print*, ""
      PRINT*, 'DYNAMIC EQUILIBRATION LOOP'
      PRINT*, 0, REAL(vacf(:,tini),sp)

      OPEN(99, file='output/vacf.dat')
      OPEN(102, FILE='output/diff.dat')
      OPEN(103, FILE='output/disp.dat')
        WRITE(99,*) 0, vacf(x,past), vacf(y,past), vacf(z,past)
        WRITE(102,*) 0, 0.0, 0.0, 0.0
        WRITE(103,*) 0, 0.0, 0.0, 0.
      CLOSE(99)
      CLOSE(102)
      CLOSE(103)
    END SUBROUTINE init_propagation


SUBROUTINE PROPAGATE(it, is_converged)

      use system, only: fluid, solid, n, node, tic
      use input, only: input_char, input_log
      use OMP_lib
      implicit none
      integer(i2b), intent(in) :: it
      real(dp) :: fermi, fractionOfParticleRemaining, scattprop, scattprop_p, n_loc(lbm%lmin:lbm%lmax)
      real(dp) :: u_star(x:z), Propagated_Quantity_loc(x:z), Free_density_loc(x:z), density_loc
      logical :: error, interfacial_loc
      logical, intent(out) :: is_converged
      integer :: temp
      integer(i4b) :: tac_loc, time_loc, count_rate, days, hours, minutes, seconds
      character(len=12) :: time_str, temp_str(4)

      error=.false.

      ! ompvar = input_char("openmpover") ! this prepares the code to parallelize over x, y or z slices depending on lb.in. TODO

!$OMP PARALLEL DO PRIVATE(i,j,k,l,ip,jp,kp,fermi,scattprop,l_inv_loc,scattprop_p,n_loc,density_loc,nature_loc,u_star) &
!$OMP PRIVATE(ip_all,jp_all,kp_all,fractionOfParticleRemaining,Propagated_Quantity_loc, Free_density_loc, interfacial_loc) &
!$OMP SHARED(nature,n,lambda,Propagated_Quantity, Free_density,l_inv,a0,c,ll,lu,lx,ly,lz,considerAdsorption,tracer) &
!$OMP SHARED(Propagated_Quantity_Adsorbed,density,Adsorbed_density,interfacial,error,alpha) &
!$OMP REDUCTION(+:vacf) &
!$OMP DEFAULT(NONE)

      do k=1,lz ! we parallelize over k. If system is 30x30x1 parallelization is useless!
        kp_all(:) = [( pbc( k+c(z,l) ,z) , l=ll, lu)]
        do j=1,ly
          jp_all(:) = [( pbc( j+c(y,l) ,y) , l=ll, lu)]
          do i=1,lx
            nature_loc = nature(i,j,k)
            if (nature_loc/=fluid) cycle
            ip_all(:) = [( pbc( i+c(x,l) ,x) , l=ll, lu)]
            u_star(:) = 0.0_dp ! the average velocity at r
            fractionOfParticleRemaining = 1.0_dp ! fraction of particles staying at r; decreases in the loop over neighbours
            n_loc(:) = n(i,j,k,:)
            density_loc = density(i,j,k)
            interfacial_loc = Interfacial(i,j,k)
            Propagated_Quantity_loc(:) = Propagated_Quantity(:,i,j,k,next)
            Free_density_loc(:) = Free_density(:,i,j,k,next)
            do l = ll+1, lu ! ll is velocity=0
              ip = ip_all(l)
              jp = jp_all(l)
              kp = kp_all(l)
              if ( nature(ip,jp,kp) /= fluid ) cycle
              fermi = 1.0_dp/(1.0_dp + calc_exp_dphi(i,j,k,ip,jp,kp))
              scattprop = calc_scattprop( n_loc(l), density_loc, a0(l), lambda, fermi) ! scattering probability at r
              fractionOfParticleRemaining = fractionOfParticleRemaining - scattprop ! what is scattered away is not found anymore at r
              u_star(:) = u_star(:) + scattprop*c(:,l)
              l_inv_loc = l_inv(l)
              scattprop_p = calc_scattprop( n(ip,jp,kp,l_inv_loc), density(ip,jp,kp), a0(l_inv_loc), lambda, 1.0_dp-fermi)
              Propagated_Quantity_loc(:) = Propagated_Quantity_loc(:) + Propagated_Quantity(:,ip,jp,kp,now)*scattprop_p
              Free_density_loc(:) = Free_density_loc(:) + Free_density(:,ip,jp,kp,now)*scattprop_p
            end do
            Propagated_Quantity(:,i,j,k,next) = Propagated_Quantity_loc(:)
            Free_density(:,i,j,k,next) = Free_density_loc(:)
            vacf(:,now) = vacf(:,now) + Propagated_Quantity(:,i,j,k,now)*u_star(:)

            ! NOW, UPDATE THE PROPAGATED QUANTITIES
            if (   (.not.interfacial_loc .and. considerAdsorption) &
              .or. (.not.considerAdsorption) )then
              Propagated_Quantity(:,i,j,k,next) = &
                Propagated_Quantity(:,i,j,k,next) + fractionOfParticleRemaining*Propagated_Quantity(:,i,j,k,now)
              Free_density(:,i,j,k,next) = &
                Free_density(:,i,j,k,next) + fractionOfParticleRemaining*Free_density(:,i,j,k,now)
            else if ( interfacial_loc .and. considerAdsorption ) then
              fractionOfParticleRemaining = fractionOfParticleRemaining&
                                            & - tracer%ka*(1 - sum(Adsorbed_density(:,i,j,k,now))*alpha) 
              Propagated_Quantity(:,i,j,k,next) = Propagated_Quantity (:,i,j,k,next) &
                + fractionOfParticleRemaining * Propagated_Quantity (:,i,j,k,now) &
                + Propagated_Quantity_Adsorbed (:,i,j,k,now) * tracer%kd
              Free_density(:,i,j,k,next) = Free_density(:,i,j,k,next) &
                + fractionOfParticleRemaining * Free_density(:,i,j,k,now) &
                + Adsorbed_density(:,i,j,k,now) * tracer%kd
              Propagated_Quantity_Adsorbed(:,i,j,k,next) = &
                Propagated_Quantity_Adsorbed(:,i,j,k,now) * (1.0_dp - tracer%kd) &
                + (1 - sum(Adsorbed_density(:,i,j,k,now))*alpha)*Propagated_Quantity(:,i,j,k,now)*tracer%ka
              Adsorbed_density(:,i,j,k,next) = &
                Adsorbed_density(:,i,j,k,now) * (1.0_dp - tracer%kd) &
                + (1 - sum(Adsorbed_density(:,i,j,k,now))*alpha)*Free_density(:,i,j,k,now)*tracer%ka
            else
              stop "OMG YOU KILLED THE QUEEN!"
            end if

            if (abs(fractionOfParticleRemaining)<epsilon(1._dp)) error=.true.
          end do
        end do
      end do

!$OMP END PARALLEL DO


      if(error) stop 'somewhere restpart is negative' ! TODO one should find a better function for ads and des, as did Benjamin for pi

      DO temp=1,3
        diff_now(temp) = diff_past(temp) + (vacf(temp,now)+vacf(temp,past))/2
        disp_now(temp) = diff_now(temp) - vacf(temp,now)*it 
        deltad(temp) = abs(disp_now(temp)-disp_past(temp))
        !deltad(temp) = abs(sum(Propagated_Quantity(temp, :,:,:,next)) - sum(Propagated_Quantity(temp,:,:,:,now)))
      END DO

      ! back to the futur: the futur is now, and reinit futur
      Propagated_Quantity(:,:,:,:,now) = Propagated_Quantity(:,:,:,:,next)
      Propagated_Quantity(:,:,:,:,next) = 0.0_dp
      Free_density(:,:,:,:,now) = Free_density(:,:,:,:,next)
      Free_density(:,:,:,:,next) = 0.0_dp
      if (considerAdsorption) then
        Propagated_Quantity_Adsorbed(:,:,:,:,now) = Propagated_Quantity_Adsorbed(:,:,:,:,next)
        Propagated_Quantity_Adsorbed(:,:,:,:,next) = 0.0_dp
        Adsorbed_density(:,:,:,:,now) = Adsorbed_density(:,:,:,:,next)
        Adsorbed_density(:,:,:,:,next) = 0.0_dp
      end if

      fa_now = sum(Adsorbed_density(:,:,:,:,now))
      deltaf = fa_now - fa_past
 
      if(modulo(it,10000)==0) then 
        CALL SYSTEM_CLOCK(tac_loc, count_rate)
        time_loc = (tac_loc-tic)/count_rate
        days = floor(time_loc/24.0/3600.0)
        time_loc = time_loc - days*24.0*3600
        hours = floor(time_loc/3600.0)
        time_loc = time_loc - hours*3600.0
        minutes = floor(time_loc/60.0)
        seconds = time_loc - minutes*60.0
        write(temp_str(1),'(I2.2)') days
        write(temp_str(2),'(I2.2)') hours
        write(temp_str(3),'(I2.2)') minutes
        write(temp_str(4),'(I2.2)') seconds
        time_str = trim(temp_str(1))//':'//trim(temp_str(2))//':'//trim(temp_str(3))//':'//trim(temp_str(4))

        OPEN(99, FILE='output/vacf.dat', ACCESS='append')
        OPEN(102, FILE='output/diff.dat', ACCESS='append')
        OPEN(103, FILE='output/disp.dat', ACCESS='append') 
            print*,it,REAL(vacf(:,now),sp), real(fa_now), time_str
            !print*, it, sum(Free_density(:,:,:,:,now)), sum(Adsorbed_density(:,:,:,:,now)),&
            !                &sum(Free_density(:,:,:,:,now)) + sum(Adsorbed_density(:,:,:,:,now))   
            WRITE(99,*) it, vacf(x,now), vacf(y,now), vacf(z,now)
            WRITE(102,*) it, diff_now(x), diff_now(y), diff_now(z)
            WRITE(103,*) it, disp_now(x), disp_now(y), disp_now(z)
        CLOSE(99)
        CLOSE(102)
        CLOSE(103)
    end if


      vacf(:,past) = vacf(:,now)
      vacf(:,now) = 0.0_dp
      diff_past(:) = diff_now(:)
      disp_past(:) = disp_now(:)
      fa_past = fa_now

    !if (it == 500 .or. it == 1000 .or. it == 2000 .or. it == 3000 .or. it == 4000&
    !        & .or. it == 5000 .or. it == 10000 .or. it == 20000) then
    !    call plot_densities_of_particles(it)
    !end if

      IF (it> tmin .and. all(deltad < delta_disp)) THEN
        print*, 'CONVERGENCE IS REACHED: delta_disp <', delta_disp
        is_converged = .true.
        print*, ""
        print*, 'disp = ', disp_now(x), disp_now(y), disp_now(z)
        !call plot_densities_of_particles(it)
    ELSE IF (it> tmin .and. all(abs(vacf)<delta_vacf)) THEN
        print*, 'CONVERGENCE IS REACHED: vacf <', delta_vacf
        print*, ""
        print*, 'disp = ', disp_now(x), disp_now(y), disp_now(z)
        is_converged = .true.
        !call plot_densities_of_particles(it)
    ELSE 
        is_converged = .false.
        if (it == tmax-tmom) then
            !call plot_densities_of_particles(it)
            print*, '!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!'
            print*, '!!!!!!!!!NO CONVERGENCE!!!!!!!!!!!'
            print*, ""
        print*, 'disp = ', disp_now(x), disp_now(y), disp_now(z)
        end if    
    END IF

    IF (ISNAN(vacf(x,past)) .or. ISNAN(vacf(y,past)) .or. ISNAN(vacf(z,past)) ) then
        is_converged = .true.
        print*, 'Moment Propagation ended: Vacf = NaN'
    ENDIF

END SUBROUTINE PROPAGATE

    ! ==============================================================================
    
    SUBROUTINE plot_densities_of_particles(it)
        integer :: i, j, k, it
        character(len=40) :: txt

        write(txt,'(I6.6)') it
        OPEN(100, FILE='output/Free_density.dat')
        OPEN(101, FILE='output/Adsorbed_density.dat')
            DO i=1,lx; DO j=1,ly; DO k=1,lz;
                WRITE(100,*)i,j,k,SUM(Free_density(:,i,j,k,now))
                WRITE(101,*)i,j,k,SUM(Adsorbed_density(:,i,j,k,now))
            END DO; END DO; END DO;
        CLOSE(100)
        CLOSE(101)
    END SUBROUTINE plot_densities_of_particles    

!==============================================================================

    SUBROUTINE plot_surf_weight
        integer :: i, j, k
        OPEN(102, FILE='output/Surf_weight.dat')
        DO i=1,lx; DO j=1,ly; DO k=1,lz;
            WRITE(102,*)i,j,k,surf_weight(i,j,k)
        END DO; END DO; END DO;
        CLOSE(102)
   END SUBROUTINE plot_surf_weight

!==============================================================================

  PURE FUNCTION CALC_EXP_DPHI( i, j, k, ip, jp, kp)
    implicit none
    REAL(DP) :: CALC_EXP_DPHI
    INTEGER(i2b), INTENT(IN) :: i, j, k, ip, jp, kp
    IF( abs(tracer%z) <= epsilon(1._dp) ) THEN
      calc_exp_dphi = 1.0_dp
    ELSE
      calc_exp_dphi = EXP( tracer%z * dphi(i,j,k,ip,jp,kp) )
    END IF
  END FUNCTION CALC_EXP_DPHI

  ! ==============================================================================

  PURE FUNCTION DPHI(i,j,k,ip,jp,kp)
    use system, only: phi, elec_slope
    implicit none
    REAL(DP) :: dphi
    integer(i2b), intent(in) :: i, j, k, ip, jp, kp
    dphi = phi(ip,jp,kp) - phi(i,j,k)
    if      (i==lx .and. ip==1) then
      dphi = dphi + elec_slope(x)*(lx+1)
    else if(i==1 .and. ip==lx) then
      dphi = dphi - elec_slope(x)*(lx+1)
    else if(j==ly .and. jp==1) then
      dphi = dphi + elec_slope(y)*(ly+1)
    else if(j==1 .and. jp==ly) then
      dphi = dphi - elec_slope(y)*(ly+1)
    else if(k==lz .and. kp==1) then
      dphi = dphi + elec_slope(z)*(lz+1)
    else if(k==1 .and. kp==lz) then
      dphi = dphi - elec_slope(z)*(lz+1)
    end if
  END FUNCTION DPHI

          ! ==============================================================================

  PURE FUNCTION CALC_LAMBDA(d)
    use system, only: kBT
    implicit none
    real(dp) :: calc_lambda
    real(dp), intent(in) :: d ! diffusion coefficient to be used to compute lambda. See review by Ladd.
    calc_lambda = 4.0_dp*d/kBT
  END FUNCTION CALC_LAMBDA

! ==============================================================================

  pure function calc_scattprop (n,rho,w,lambda,fermi)
    implicit none
    real(dp) :: calc_scattprop
    real(dp), intent(in) :: n, rho, w, lambda, fermi
    calc_scattprop = n/rho - w + lambda*w*fermi
  end function calc_scattprop

! ==============================================================================

    SUBROUTINE DEALLOCATE_PROPAGATED_QUANTITY
      if (allocated(Propagated_Quantity)) deallocate(Propagated_Quantity)
      if (allocated(Propagated_Quantity_Adsorbed)) deallocate(Propagated_Quantity_Adsorbed)
    END SUBROUTINE DEALLOCATE_PROPAGATED_QUANTITY

! ==============================================================================

    SUBROUTINE test_and_allocate_what_is_needed_for_moment_propagation
      implicit none
      ! Propagated_Quantity is the probability vector of arriving at r at time t
      IF(ALLOCATED(Propagated_Quantity)) STOP 'Propagated quantity should not be allocated in init_moment_propagation'
      ALLOCATE(Propagated_Quantity(x:z,lx,ly,lz,now:next), source=0.0_dp)
      IF(ALLOCATED(Free_density)) STOP 'Free density should not be allocated in init_moment_propagation'
      ALLOCATE(Free_density(x:z,lx,ly,lz,now:next), source=0.0_dp)
      
      !if (considerAdsorption) then
        IF(ALLOCATED(Propagated_Quantity_Adsorbed)) STOP 'Propagated_Quantity_Adsorbed should not be allocated yet'
        ALLOCATE(Propagated_Quantity_Adsorbed(x:z,lx,ly,lz,now:next), source=0.0_dp)
        IF(ALLOCATED(Adsorbed_density)) STOP 'Adsorbed density should not be allocated yet'
        ALLOCATE(Adsorbed_density(x:z,lx,ly,lz,now:next), source=0.0_dp)
      !end if
    END SUBROUTINE test_and_allocate_what_is_needed_for_moment_propagation

END MODULE MOMENT_PROPAGATION
