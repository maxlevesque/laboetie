! Here we initiate the charged densities (solutes + and -)

subroutine charges_init

    use precision_kinds, only: i2b, dp
    use system, only: lambda_D, c_plus, c_minus, phi, charge_distrib, sigma, node,&
                       fluid, solid, anormf0, bjl, kBT, rho_0, D_plus, D_minus, supercell
    use constants, only: pi, x, y, z
    use module_input, only: getinput
    use myallocations

    implicit none
    real(dp), allocatable, dimension(:,:,:) :: c_plusT
    integer  :: count_solid, count_fluid, count_solid_int, geometrie, i,j,k, lx, ly, lz, ENNE,&
                capacitor
    real(dp) :: sigma_solid, sigma_fluid, SF, PrefactLP1, PrefactLP2, PrefactLP3, Alpha
    real(dp) :: in_c_plus_solid, in_c_plus_fluid, in_c_minus_solid, in_c_minus_fluid
    REAL(dp), PARAMETER :: eps=EPSILON(1._dp)

     open(279, file='output/c_plus_alongZCHARGESINIT.dat')
     open(277, file='output/c_at_wall1.dat')
     open(276, file='output/c_at_wall2.dat')


    !
    ! Is there any charge into the solute ?
    !
    bjl = getinput%dp('bjl',0._dp)
    capacitor = getinput%int('capacitor',0) ! 1 = true 0 = false

    IF( bjl <= eps ) RETURN

    lx = supercell%geometry%dimensions%indiceMax(x)
    ly = supercell%geometry%dimensions%indiceMax(y)
    lz = supercell%geometry%dimensions%indiceMax(z)

    geometrie = getinput%int('geometryLabel',-1) ! Ade: 15/03/2017
    sigma = getinput%dp('sigma',0._dp)
    lambda_d = getinput%dp('lambda_D',0._dp)

    if( abs(lambda_D) <= epsilon(1._dp) ) then
      print*, 'salt free fluid'
      rho_0 = 0.0_dp
    else
      rho_0 = 1.0_dp / (4.0_dp*pi*bjl*lambda_D**2)
    end if

    ! count of solid, fluid and interfacial nodes
    count_solid = count(node%nature==solid)
    count_solid_int = count( node%nature == solid .and. node%isInterfacial )
    count_fluid = count(node%nature==fluid)

    ! read where are distributed the charges
    ! call read_charge_distrib

    !---------------------------Ade----------------------------------
    ! We moved the allocation of the following table outside the if statement underneath
    !if (.not. allocated(c_plus)) call allocateReal3D(c_plus)
    !if (.not. allocated(c_minus)) call allocateReal3D( c_minus)
    !---------------------------Ade----------------------------------
    if (.not. allocated(c_plusT)) call allocateReal3D(c_plusT)

    IF( ABS(sigma) > EPSILON(1._dp) ) THEN ! Ade : this if statement should be removed 27/01/2017

        ! init ion (solute) concentrations
        if (.not. allocated(c_plus)) call allocateReal3D(c_plus)
        if (.not. allocated(c_minus)) call allocateReal3D( c_minus)

        ! init potential
        call allocateReal3D( phi) !allocate( phi(lx,ly,lz), source=0.0_dp )
        phi = 0._dp
        charge_distrib = getinput%char("charge_distrib")
        if( charge_distrib(1:3) /= 'int' .and. charge_distrib(1:3)/='sol') stop 'charge_distrib can only be int or sol for now'

        ! distribute charge, depending on where user asked
        if( charge_distrib(1:3) == 'sol') then
          if(count_solid/=0) then
            sigma_solid = sigma / count_solid ! charges distributed in all solid nodes
          else
            sigma_solid = 0
          end if

        else if( charge_distrib(1:3) == 'int') then
          if(count_solid_int/=0) then
            sigma_solid = sigma / count_solid_int ! charges distributed in interfacial solid nodes only
          else
            sigma_solid = 0
          end if
        end if

        if(count_fluid/=0) then
          sigma_fluid = -1.0_dp * sigma / count_fluid ! charges distributed in all fluid nodes
        else
          sigma_fluid = 0
        end if

        in_c_plus_solid  = +0.5_dp*sigma_solid;
        in_c_minus_solid = -0.5_dp*sigma_solid;

        if( sigma_fluid > 0.0_dp ) then
              in_c_plus_fluid  = 0.5_dp*rho_0 + sigma_fluid
              in_c_minus_fluid = 0.5_dp*rho_0
        else
          in_c_plus_fluid  = 0.5_dp*rho_0
          in_c_minus_fluid = 0.5_dp*rho_0 - sigma_fluid
        end if

        ! Ade : modification 15/03/2017
       if(capacitor.NE.1) then
          if( geometrie == 1 ) then   ! slit pore geometry
              Alpha = getinput%dp('Alpha',0._dp)  ! Attention!!!!!! This is dangerous. We should probably  
                                                  ! compute its value in another subroutine
              ENNE = lz-2 ! Nb of fluid nodes
              if(mod(lz,2) == 0) then
                 SF = real(lz)/2
              else
                 SF = real(lz)/2 + 0.5
              endif
              if( lambda_d > EPSILON(1._dp)) then ! Low potential condition - salt added
                  !SurfArea = (lx*ly)_dp
                  PrefactLP1 = 1/(8*PI*bjl*lambda_d**2) 
                  PrefactLP2 = 4*PI*sigma*bjl*lambda_d/(2*lx*ly)
                  PrefactLP3 = sinh(ENNE / ( 2*lambda_d) )
                  do i = 1, lx    
                    do j = 1, ly   
                      do k = 2, lz-1 ! first and last node are solid
                        c_plusT(i,j,k) = PrefactLP1 * ( 1 - PrefactLP2 * cosh( (k-SF)/lambda_d )/(PrefactLP3) ) !* SurfArea
                      end do
                    end do
                  end do
              else      ! Normal c_plus density - no salt
                  do i = 1, lx      
                    do j = 1, ly  
                      do k = 2, lz-1 ! first and last node are solid
                         c_plusT(i,j,k) = (Alpha**2/( 2*PI*bjl * ( cos(Alpha*(k-SF)) )**2 )) !* SurfArea
                      end do
                    end do
                  end do
              endif
          endif
        end if

        ! Ade : end modification 15/03/2017
        
        if(abs(sigma)>epsilon(1._dp)) then
          if(charge_distrib(1:3)=='int') then
            print*,'The total charge is set ONLY on the solid nodes at the interface (',count_solid_int,'/',count_solid,')'
            print*,'Internal sites (at the interface) = ',count_solid_int,', charge per link =',sigma_solid
            print*,'External sites = ',count_fluid,' charge per link =',sigma_fluid
          else if(charge_distrib(1:3)=='sol') then
            print*,'Internal sites =',count_solid,' charge per link =',sigma_solid
            print*,'External sites =',count_fluid,' charge per link =',sigma_fluid
            stop 'only surface charge is implemented for now in charge_init.f90'
          else
            stop 'pb in charges_init.f90'
          end if

          print*,'Salt concentration ',0.5*rho_0
          print*,'Init density values :'
          print*,'p_solid =',in_c_plus_solid
          print*,'p_fluid =',in_c_plus_fluid
          print*,'m_solid =',in_c_minus_solid
          print*,'m_fluid =',in_c_minus_fluid
          print*,'*********************************************************************'

          print*,'ATTENTION ONLY SURFACE CHARGE IS OK FOR NOW'
        end if

        where(node%nature==solid .and. node%isInterfacial )
          c_plus = in_c_plus_solid
          c_minus = in_c_minus_solid
        else where(node%nature==solid .and. .not. node%isInterfacial )
          c_plus = 0.0_dp
          c_minus = 0.0_dp
        else where(node%nature==fluid)
            c_plus = in_c_plus_fluid
            c_minus = in_c_minus_fluid
        end where

        ! Ade : modification 24/03/2017
        ! This modication was done in order to simulate a capacitor. A wall will have then 
        ! negative charges, whilst the other will have positive charges
        if (capacitor.EQ.1) then
            do j=1,ly
              c_plus(:,:,1) = in_c_plus_solid/2
              c_minus(:,:,1) = in_c_plus_solid/2
            print*, ' ------------------------------------ '
            print*, 'c_plus wall1 = ', c_plus
              c_plus(:,:,lz) = -in_c_plus_solid/2
              c_minus(:,:,lz) = -in_c_plus_solid/2
            print*, ' ------------------------------------ '
            print*, 'c_plus wall2 = ', c_plus
            end do
        !    do k = 1,lz
        !      if(k==1) then
        !        where(node%nature==solid .and. node%isInterfacial )
        !          c_plus = in_c_plus_solid/2    ! negative c_plus
        !          c_minus = in_c_minus_solid/2
        !        else where(node%nature==solid .and. .not. node%isInterfacial )
        !          c_plus = 0.0_dp
        !          c_minus = 0.0_dp
        !        end where
        !          print*, ' ------------------------------------ '
        !          print*, 'c_plus wall1 = ', c_plus
        !      else if (k==lz) then
        !        where(node%nature==solid .and. node%isInterfacial )
        !          c_plus = -in_c_plus_solid/2 ! positive c_plus
        !          c_minus = -in_c_minus_solid/2
        !        else where(node%nature==solid .and. .not. node%isInterfacial )
        !          c_plus = 0.0_dp
        !          c_minus = 0.0_dp
        !        end where
        !          print*, ' ------------------------------------ '
        !          print*, 'c_plus wall2 = ', c_plus
        !      endif
        !    end do
        end if 
        ! Ade : modification 24/03/2017

        ! Ade : modification 15/03/2017
        ! Attention!!! you should write also the corresponding c_minus
        if (capacitor==1) then
            where(node%nature==fluid)
              c_plus = 0.0_dp
              c_minus = 0.0_dp
            end where
        else if (geometrie == 1 .and. capacitor.NE.1) then
              where(node%nature==fluid)
                c_plus = c_plusT
              end where
        end if
        ! Ade : end modification 15/03/2017
        
        
        anormf0 = 4.0_dp*pi*bjl*kBT/2.0_dp *sum(abs(c_plus)+abs(c_minus))

        ! ---------------------- Ade ---------------------------
        !print*, 2, SUM(c_plus(:,:,2))
        ! Ade : this is a debugging test
        ! ---------------------- Ade ---------------------------

        DO k=supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
           WRITE(279,*) k, SUM(c_plus(:,:,k))
           if(k==1) then
              write(277,*) k, SUM(c_plus(:,:,k))
           else if (k==lz) then
              write(276,*) k, SUM(c_plus(:,:,k))
           end if 
        ENDDO
        close(279)

        ! TODO call charge_test

        ! read diffusion coefficients of solute + and solute -
        d_plus = getinput%dp('D_plus',0._dp)
        if( D_plus < 0.0_dp ) stop 'D_plus <0. critical.'
        d_minus = getinput%dp('D_minus',0._dp)
        if( D_minus < 0.0_dp ) stop 'D_minus <0. critical.'
    END IF
end subroutine charges_init
