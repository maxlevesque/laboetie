! Here is a smolu with diffusion + electrostatic part only
! The objective is to end up in a situation where equilibrium
! distribution of solutes (c_plus and c_minus) is found,
! ie in which the total flux of solutes is zero.

subroutine just_eq_smolu

    use precision_kinds, only: dp
    use system, only: D_plus, D_minus, node,&
            D_equil, time, fluid, solid, phi, c_plus, c_minus,&
            rho_0, sigma, pbc, supercell, x, y, z
    use mod_lbmodel, only: lbm
    use myAllocations, only: allocateReal3D
    
    implicit none
    
    integer :: iter, max_iter
    integer :: fluidNodesCount ! count of fluid nodes in the supercell
    real(dp), allocatable :: flux_site_plus(:,:,:), flux_site_minus(:,:,:)
    integer :: i, j, k, ip, jp, kp, l, lx, ly, lz ! dummy
    real(dp) :: exp_dphi, exp_min_dphi ! Exp[phi site1 - phi site 2] and 1/Exp
    real(dp) :: flux_link_plus, flux_link_minus
    real(dp) :: tot_diff_plus, tot_diff_minus ! total flux of + and - solutes init at high values
    real(dp) :: eD_plus, eD_minus ! effective D_plus and D_minus
    real(dp), parameter :: convergence_criteria = 5.e-6 ! magic number !?
    real(dp) :: FactMinus, FactPlus ! dummy variables

    call allocateReal3D( flux_site_minus )
    call allocateReal3D( flux_site_plus  )

    fluidNodesCount = count(node%nature==fluid)

    !
    ! init
    !
    
    ! diffusion coefficients of solutes are read in input file.
    ! if we're in the first steps of equilibration, it speeds up convergence
    ! to use a smaller diffusion coefficient.
    ! i think there is a bug in C code here,
    ! especially considering the high amount of magic numbers and +- convention for time.
    tot_diff_plus = convergence_criteria
    tot_diff_minus = convergence_criteria


    if( time > 0 .or. time < -D_equil ) then
        error stop 'in just_eq_smolu.f90. Should be accessed in equilibration steps only'
    else if( time < -D_equil+100) then ! 100 first steps
        eD_plus = 0.1_dp*D_plus
        eD_minus = 0.1_dp*D_minus
        max_iter = 10
    else if( time < -D_equil +500) then
        eD_plus = D_plus
        eD_minus = D_minus
        max_iter = 1
    else

        if( D_plus < 0.03_dp) then
            eD_plus = 0.03_dp ! too a small D is not effective for this step. no physic is associated to the "dynamic" of this step.
        else
            eD_plus = D_plus
        end if

        if( D_minus < 0.03_dp) then
            eD_minus = 0.03_dp
        else
            eD_minus = D_minus
        end if
        max_iter = 1

    end if

    print*,'D_plus, D_minus, max_iter ', eD_plus, eD_minus, max_iter

    lx = supercell%geometry%dimensions%indiceMax(x)
    ly = supercell%geometry%dimensions%indiceMax(y)
    lz = supercell%geometry%dimensions%indiceMax(z)

    !
    ! init convergence iterations
    !
    iter = 0

    do while( tot_diff_minus + tot_diff_plus > convergence_criteria .and. iter < max_iter )
        ! tot_diff_minus > convergence_criteria .and. tot_diff_plus > convergence_criteria would make more sense to me.

        ! flux_site_plus or minus is what makes the system evoluate.
        flux_site_plus = 0.0_dp
        flux_site_minus = 0.0_dp

        ! for all sites
        do k = 1, lz
            do j = 1, ly
                do i = 1, lx

                    ! the flux from site (i,j,k) is 0 if (i,j,k) is a solid node
                    if( node(i,j,k)%nature == solid ) cycle

                    ! for each neighbouring site, i.e., accessible site
                    do l = lbm%lmin +1, lbm%lmax, 2
                    
                        ! periodic boundary conditions to neighbours
                        ip = pbc( i+ lbm%vel(l)%coo(x) ,x)
                        jp = pbc( j+ lbm%vel(l)%coo(y) ,y)
                        kp = pbc( k+ lbm%vel(l)%coo(z) ,z)
                        if( node(ip,jp,kp)%nature == solid) cycle ! no flux of matter toward solid nodes of course

                        ! compute the difference in potential between site (i,j,k) and site (ip,jp,kp)
                        exp_dphi = exp( phi(ip,jp,kp) - phi(i,j,k) ) ! be carefull to sign
                        exp_min_dphi = 1.0_dp / exp_dphi

                        ! flux due to electrostatic and density gradients inside link i,j,k <-> ip,jp,kp
                        flux_link_plus  = 0.5_dp * (1.0_dp+ exp_min_dphi)*( c_plus(ip,jp,kp) * exp_dphi     - c_plus (i,j,k) )
                        flux_link_plus = flux_link_plus * eD_plus / lbm%vel(l)%delta
                        flux_link_minus = 0.5_dp * (1.0_dp+ exp_dphi    )*( c_minus(ip,jp,kp) * exp_min_dphi - c_minus(i,j,k) )
                        flux_link_minus = flux_link_minus * eD_minus / lbm%vel(l)%delta

                        ! update flux on each site accordingly to flux inside link involving site
                        flux_site_plus(i,j,k) = flux_site_plus(i,j,k) + flux_link_plus
                        flux_site_minus(i,j,k) = flux_site_minus(i,j,k) + flux_link_minus
                        flux_site_plus(ip,jp,kp) = flux_site_plus(ip,jp,kp) - flux_link_plus
                        flux_site_minus(ip,jp,kp) = flux_site_minus(ip,jp,kp) - flux_link_minus
                    
                    end do ! l

                end do ! i
            end do ! j
        end do ! k

        ! Check that no matter has been lost. It should be transfered only. Thus, global concentrations must not change, ie sum of flux over all sites should be 0.
        if( abs(sum(flux_site_plus)) > 1.e-12 .or. abs(sum(flux_site_minus)) > 1.e-12) then
            print*,'abs(sum(flux_site_plus))=',abs(sum(flux_site_plus))
            print*,'abs(sum(flux_site_minus))=',abs(sum(flux_site_minus))
            error stop 'the sum of all flux does not add up. problem in just_eq_smolu.f90'
        end if

        !
        ! update local concentrations (smolushowski part)
        !
        where( node%nature == fluid ) ! select all nodes (i,j,k) that are fluid
            c_plus = c_plus + flux_site_plus
            c_minus = c_minus + flux_site_minus
        end where

        !
        ! compute the total flux in this equilibration step one wants to minimize.
        !
        if( sigma/=0 .and. eD_plus/=0.0_dp .and. eD_minus/=0.0_dp ) then
            ! the sum of all flux
            FactPlus  = sum(c_plus)  / fluidNodesCount + rho_0
            FactMinus = sum(c_minus) / fluidNodesCount + rho_0
            tot_diff_plus  = sqrt(sum(flux_site_plus**2,mask=(node%nature==fluid))) &
                / fluidNodesCount / (0.5_dp*FactPlus*eD_plus ) / sigma
            tot_diff_minus = sqrt(sum(flux_site_minus**2,mask=(node%nature==fluid)))&
                / fluidNodesCount / (0.5_dp*FactMinus*eD_minus) / sigma
        end if

        ! increment iteration count
        iter = iter +1

    end do ! while loop about convergence on tot_diff_minus+tot_diff_plus

    deallocate(flux_site_minus)
    deallocate(flux_site_plus)

end subroutine just_eq_smolu
