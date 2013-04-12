! Here we evolve (propagate) the population n accordingly to the velocities l
subroutine propagation
  use precision_kinds, only: i2b, dp
  use system, only: n, pbc, supercell
  use constants, only: x, y, z
  use mod_lbmodel, only: lbm
  implicit none
  integer(i2b) :: i, j, k, l, ip, jp, kp
  real(dp), allocatable, dimension(:,:,:) :: old_n
  integer(i2b) :: lx, ly, lz
  lx = supercell%geometry%dimensions%indiceMax(x)
  ly = supercell%geometry%dimensions%indiceMax(y)
  lz = supercell%geometry%dimensions%indiceMax(z)
  allocate(old_n(lx,ly,lz))
  call boundPM
  ! for each velocity, find the futur node and put it population at current time
  do l= lbm%lmin, lbm%lmax
    old_n = n(:,:,:,l) ! backup population before propagation
    ! propagate for each starting node i,j,k to ip, jp, kp
    do k=1,lz
      kp=pbc(k+lbm%vel(l)%coo(z),z)
      do j=1,ly
        jp=pbc(j+lbm%vel(l)%coo(y),y)
        do i=1,lx
          ip=pbc(i+lbm%vel(l)%coo(x),x)
          n (ip,jp,kp,l) = old_n (i,j,k) ! evolve the system (propagate) ie evolve population of (i,j,k) to (ip,jp,kp)
        end do
      end do
    end do
  end do
end subroutine propagation





!/***************************************************************/
!/**   This subroutine exchange the velocities of two*/
!/**   neigbouring sites if one is a FLUID site an the*/
!/**   other is a SOLID state.*/
!/**   It prepares the 'n's to be propagated all in the SAME way*/
!/**/*/
!As the lattice Boltzmann method is a kinetic method, macroscopic boundary con-
!ditions do not have direct equivalents. They have to be replaced by appropriate mi-
!croscopic rules which induce the desired macroscopic behavior. The easiest solution
!for introducing solid walls (i.e. a no-slip boundary condition) is the introduction of
!the bounce back rule on wall nodes
!   f_out (x,l,t) = f_in(x,-l,t), with x+l ∈ wall,
!This rule rotates the distribution functions on the wall node and
!thus they return back to the fluid with opposite momentum in the next time step. This
!results in zero velocities at the wall (which is located half-way between the last fluid
!cell and the first wall node) and ensures that there is no flux across the wall. This is
!equivalent to the macroscopic no-slip boundary condition.

SUBROUTINE BOUNDPM
  use precision_kinds
  use system, only: pbc, supercell, n,solid,fluid
  use constants, only: x, y, z
  use mod_lbmodel, only: lbm
  implicit none
  integer(i2b) :: i, j, k, l, w, ip, jp, kp
  real(dp) :: tmp
  do i= 1, supercell%geometry%dimensions%indiceMax(x)
    do j= 1, supercell%geometry%dimensions%indiceMax(y)
      do k= 1, supercell%geometry%dimensions%indiceMax(z)
        do l= lbm%lmin, lbm%lmax, 2
          ip= pbc( i+ lbm%vel(l)%coo(x) ,x)
          jp= pbc( j+ lbm%vel(l)%coo(y) ,y)
          kp= pbc( k+ lbm%vel(l)%coo(z) ,z)
          if( supercell%node(i,j,k)%nature /= supercell%node(ip,jp,kp)%nature ) then
            w = lbm%vel(l)%inv
            tmp = n(i,j,k,l)
            n(i,j,k,l) = n(ip,jp,kp,w)
            n(ip,jp,kp,w) = tmp
          end if
        end do
      end do
    end do
  end do
END SUBROUTINE BOUNDPM
