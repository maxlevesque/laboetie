! Here we evolve (propagate) the population n accordingly to the velocities l

SUBROUTINE PROPAGATION

  use precision_kinds, only: i2b, dp
  use system, only: n, lx, ly, lz, c, NbVel, plusx, plusy, plusz
  use constants, only: x, y, z

  implicit none
  integer(i2b) :: i, j, k, l, ip, jp, kp
  real(dp), dimension(lx,ly,lz) :: old_n

  call boundPM

  ! for each velocity, find the futur node and put it population at current time
  do l= 1, NbVel
    old_n = n(:,:,:,l) ! backup population before propagation
    ! propagate for each starting node i,j,k to ip, jp, kp
    do k=1,lz
      kp=plusz(k+c(z,l))
      do j=1,ly
        jp=plusy(j+c(y,l))
        do i=1,lx
          ip=plusx(i+c(x,l))
          n (ip,jp,kp,l) = old_n (i,j,k) ! evolve the system (propagate) ie evolve population of (i,j,k) to (ip,jp,kp)
        end do
      end do
    end do
  end do

END SUBROUTINE PROPAGATION





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
  use system, only: lx, ly, lz, nbvel, c, plusx, plusy, plusz, inside, vel_inv, n,solid,fluid
  use constants, only: x, y, z
  implicit none
  integer(i2b) :: i, j, k, l, w, ip, jp, kp
  real(dp) :: tmp
  do i= 1, lx
    do j= 1, ly
      do k= 1, lz
        do l= 1, nbvel, 2
          ip= plusx( i+ c(x, l))
          jp= plusy( j+ c(y, l))
          kp= plusz( k+ c(z, l))
          if( inside(i,j,k) /= inside(ip,jp,kp) ) then
            w = vel_inv(l)
            tmp = n(i,j,k,l)
            n(i,j,k,l) = n(ip,jp,kp,w)
            n(ip,jp,kp,w) = tmp
          end if
        end do
      end do
    end do
  end do
END SUBROUTINE BOUNDPM
