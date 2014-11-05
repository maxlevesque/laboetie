! Here we compute the density rho(i,j,k) at each node i,j,k.
! It is the sum over all velocity component l at each node.

subroutine comp_rho

  use system, only:  n, supercell, node
  implicit none

  ! new density is the sum over all velocities of the population at each node
!  rho = sum( n, 4)
  node%solventDensity = sum(n,4)
  !  print*,'Total density calculated in comp_rho.f90 = ', sum(rho)/(lx*ly*lz)

end subroutine comp_rho
