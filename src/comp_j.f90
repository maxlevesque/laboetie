! Here we compute the momentum density j. More precisely its three components against x, y and z.

subroutine comp_j

  use precision_kinds, only: i2b, dp
  use system, only: n, supercell !, f_ext, solid
  ! n is the population in direction 1, 2 and 3 of velocity discrete value l (1:19)
  ! c is the quadrature in l mathematica direction
  ! j is the momentum density in direction x, ie jx, is   jx(i,j,k) = sum_l c_x(l) * n(i,j,k,l)    where c_x(l) = c(1,l)
  use constants, only: x, y, z
  use mod_lbmodel, only: lbm

  implicit none
  integer(kind=i2b) :: i, j, k ! dummy
  do concurrent( k= supercell%geometry%dimensions%indiceMin(z):supercell%geometry%dimensions%indiceMax(z),&
                 j= supercell%geometry%dimensions%indiceMin(y):supercell%geometry%dimensions%indiceMax(y),&
                 i= supercell%geometry%dimensions%indiceMin(x):supercell%geometry%dimensions%indiceMax(x) )

! Conventional definition :
!        jx(i,j,k) = sum( n(i,j,k,:) * lbm%vel(:)%coo(x) )
!        jy(i,j,k) = sum( n(i,j,k,:) * lbm%vel(:)%coo(y) ) 
!        jz(i,j,k) = sum( n(i,j,k,:) * lbm%vel(:)%coo(z) ) 

! The historical definition (Capuani & Frenkel's code).
! Much better than the naive conventional definition in terms of numerical stability.
        supercell%node(i,j,k)%solventFlux(x) = 0.5_dp*(supercell%node(i,j,k)%solventFlux(x) + sum(n(i,j,k,:)*lbm%vel(:)%coo(x)))
        supercell%node(i,j,k)%solventFlux(y) = 0.5_dp*(supercell%node(i,j,k)%solventFlux(y) + sum(n(i,j,k,:)*lbm%vel(:)%coo(y)))
        supercell%node(i,j,k)%solventFlux(z) = 0.5_dp*(supercell%node(i,j,k)%solventFlux(z) + sum(n(i,j,k,:)*lbm%vel(:)%coo(z)))

! see Ladd and Verberg, J. Stat. Phys. 104, 1191 (2001) at page 1211 Eq. 29 :
! The forcing term is expanded in a power series in the particle velocity. The zeroth
! and first moments are given by conservation. The second moment is usually neglected,
! but in that case the external force contributes indue terms to the momentum flux.
! One may redefine the momentum density, and that's what we do here
! TODO UNSURE VERIFY THE FOLLOWING ! OCTOBER 2012
!        jx(i,j,k) = sum( n(i,j,k,:) * lbm%vel(:)%coo(x) ) + f_ext(x)/2.0_dp
!        jy(i,j,k) = sum( n(i,j,k,:) * lbm%vel(:)%coo(y) ) + f_ext(y)/2.0_dp
!        jz(i,j,k) = sum( n(i,j,k,:) * lbm%vel(:)%coo(z) ) + f_ext(z)/2.0_dp

  end do

! USELESS FOR NOW THAT BOUNDPM HAS BEEN CHECKED CAREFULLY
!  where( supercell%node%nature == solid )
!    jx = 0.0_dp
!    jy = 0.0_dp
!    jz = 0.0_dp
!  end where

end subroutine comp_j
