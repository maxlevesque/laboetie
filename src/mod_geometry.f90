MODULE GEOMETRY
  use precision_kinds, only: dp, i2b
  use system, only: lx, ly, lz, fluid, solid, inside
  implicit none
  private
  public construct_slit, construct_cylinder, construct_cc, construct_disc_benichou

CONTAINS




SUBROUTINE CONSTRUCT_SLIT
  integer(i2b) :: mi, ma
  mi = lbound(inside,3)
  ma = ubound(inside,3)
  inside = fluid
  inside( :, :, mi) = solid ! the lower bound of the thrid dimension of inside is solid
  inside( :, :, ma) = solid ! so is the upper bound
END SUBROUTINE CONSTRUCT_SLIT







SUBROUTINE CONSTRUCT_CC
  implicit none
  integer(i2b) :: i,j,k ! dummy
  if( lx /= ly .or. lx /= lz ) stop 'with wall = 3, i.e. cfc cell, the supercell should be cubic with lx=ly=lz'
  inside = fluid
  do concurrent( i=1:lx, j=1:ly, k=1:lz )
    if( is_in_solid_sphere(i,j,k) ) then
      inside(i,j,k) = solid
    else
      inside(i,j,k) = fluid
    end if
  end do
  contains
    PURE LOGICAL FUNCTION IS_IN_SOLID_SPHERE(i,j,k)
      implicit none
      integer(i2b), intent(in) :: i,j,k
      real(dp), dimension(9) :: distances
      distances(1) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/0.0_dp, 0.0_dp, 0.0_dp/) )
      distances(2) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp), 0.0_dp, 0.0_dp/) )
      distances(3) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/0.0_dp, real(ly,dp), 0.0_dp/) )
      distances(4) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/0.0_dp, 0.0_dp, real(lz,dp)/) )
      distances(5) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp), real(ly,dp), 0.0_dp/) )
      distances(6) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp), 0.0_dp, real(lz,dp)/) )
      distances(7) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/0.0_dp, real(ly,dp), real(lz,dp)/) )
      distances(8) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp), real(ly,dp), real(lz,dp)/) )
      distances(9) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp)/2.0_dp, real(ly,dp)/2.0_dp, real(lz,dp)/2.0_dp/) ) ! center
      ! put distances to all of corners of the cube and to its center (=9 points) in an array
      if ( any(distances <= real(lx,dp)*sqrt(3.0_dp)/4.0_dp) ) then
        is_in_solid_sphere = .true.
      else
        is_in_solid_sphere = .false.
      end if
    END FUNCTION IS_IN_SOLID_SPHERE
END SUBROUTINE CONSTRUCT_CC







SUBROUTINE CONSTRUCT_CYLINDER
  implicit none
  real(dp) :: radius ! radius of cylinder
  real(dp), dimension(2) :: rnode, rorigin ! coordinates of each node and center of cylinder in x,y coordinates
  integer(i2b) :: i, j ! dummy

  if( lx /= ly) stop 'wall=2 is for cylinders, which should have same lx and ly'
  rorigin = [ real(lx+1,dp)/2.0_dp, real(ly+1,dp)/2.0_dp ]
  radius = real(lx-1,dp)/2.0_dp

  do i = 1, lx
    do j = 1, ly
      rnode = [real(i,dp),real(j,dp)] - rorigin
      if( norm2(rnode) >= radius ) then ! = radius is important because without it one has exists
        inside(i,j,:) = solid
      else
        inside(i,j,:) = fluid
      end if
    end do
  end do
END SUBROUTINE CONSTRUCT_CYLINDER






SUBROUTINE CONSTRUCT_DISC_BENICHOU
! this program computes the time-dependent diffusion coefficient
! of a two-dimensional system. it is a test of Olivier Benichou's
! problem as expressed during the unformal discussion in PECSA
! to present him the numerical results of LB with sorption.
! The system is a circle, in which four entrances (exit) are found
! at 0, 3, 6 and 9".
!
!ooooo ooooo
!ooo     ooo
!oo       oo
!o         o
!o         o
!           
!o         o
!o         o
!oo       oo
!ooo     ooo
!ooooo ooooo

! Note the four exits at 0, 3, 6 and 9".
  real(dp) :: radius ! radius of cylinder
  real(dp), dimension(2) :: rnode, rorigin ! coordinates of each node and center of cylinder in x,y coordinates
  integer(i2b) :: i, j ! dummy

  if( lx /= ly) stop 'wall=2 is for cylinders, which should have same lx and ly'
  if( mod(lx,2) == 0 ) stop 'lx should be odd'
  rorigin = [ real(lx+1,dp)/2.0_dp, real(ly+1,dp)/2.0_dp ]
  radius = norm2( [1,(lx+1)/2] - rorigin) ! take great care as this is lx/2 only if rorigin falls in a lattice point.
  do i = 1, lx
    do j = 1, ly
      rnode = [real(i,dp),real(j,dp)] - rorigin
      if( norm2(rnode) > radius ) then
        inside(i,j,:) = solid
      else
        inside(i,j,:) = fluid
      end if
    end do
  end do

END SUBROUTINE CONSTRUCT_DISC_BENICHOU




END MODULE GEOMETRY
