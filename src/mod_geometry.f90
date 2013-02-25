MODULE GEOMETRY
  use precision_kinds ! all precision kinds defined in the dedicated module
  use system, only: lx, ly, lz, fluid, solid, inside
  implicit none
!  private
!  public construct_slit, construct_cylinder, construct_cc, construct_disc_benichou

CONTAINS



SUBROUTINE CONSTRUCT_TUBE_WITH_VARYING_DIAMETER
! construct the system described in Makhnovskii et al., Chem. Phys. 367, 110 (2010)
! 4 effective variables: length of big pore lw, narrow pore ln, radius of each, R and a, respectively.
! the tube is along x, with infinity sym along z (ie lz = 1 ok)
  real(dp) :: R, a, ln, lw, lw_o_lx, a_o_R
  real(dp), dimension(2) :: r0, rjk
  integer(i2b) :: i, j, k
  logical :: is_here
  if(ly/=lz) stop 'ly should be equal to lz cause it is a tube'
  inside = solid
  inquire(file='tube.in',exist=is_here)
  if(.not.is_here) stop 'tube.in cannot be read'
  open(unit=12,file='tube.in')
  read(12,*) lw_o_lx, a_o_R
  R = real(ly-1)/2.0
  a = a_o_R*R
  r0 = [real(ly+1)/2.,real(lz+1)/2.]
  do concurrent(i=1:lx, j=1:ly, k=1:lz)
    rjk = [real(j),real(k)]
    if(i>lw_o_lx*lx) then
      if( norm2( rjk-r0 ) < a ) then
        inside(i,j,k) = fluid
      else
        inside(i,j,k) = solid
      end if
    else if(i<=lw_o_lx*lx) then
      if( norm2( rjk-r0 ) < R ) then
        inside(i,j,k) = fluid
      else
        inside(i,j,k) = solid
      end if
    end if
  end do
END SUBROUTINE CONSTRUCT_TUBE_WITH_VARYING_DIAMETER






SUBROUTINE CONSTRUCT_PLANES_WITH_VARIOUS_RADIUS_2D
! construct the system described in Makhnovskii et al., Chem. Phys. 367, 110 (2010)
! 4 effective variables: length of big pore lw, narrow pore ln, radius of each, R and a, respectively.
! the tube is along x, with infinity sym along z (ie lz = 1 ok)
  real(dp) :: R, a, ln, lw
  integer(i2b) :: i
!  R = real(ly-2,dp)/2._dp
!  a = R/2._dp
!  lw = real(lx,dp)/2._dp
!  ln = lw
  inside = fluid
  inside(:,1,:) = solid
  inside(:,ly,:) = solid
  do i=1,lx
    if( i>lx/2 ) then
      inside(i,2:1+(ly-2)/4,:) = solid
      inside(i,ly-(ly-2)/4:ly-1,:) = solid
    end if
  end do
END SUBROUTINE CONSTRUCT_PLANES_WITH_VARIOUS_RADIUS_2D


SUBROUTINE CONSTRUCT_SINUSOIDAL_WALLS_2D
  ! tunnel along lx, ly is the width of the tunnel
  real(dp) :: x, y
  real(dp) :: a, b
  real(dp), parameter :: pi = acos(-1.0_dp)
  integer(i2b) :: i, j, ninty, mirror, midly
  ! ly should be odd, so that the middle of ly is on a node (w(x)=0)
  if( mod(ly,2)==0 ) stop 'ly should be odd'
  midly = (ly+1)/2
  b = 2./3.*midly
  a = midly - b - 1
  if( a<=0 .or. b<=0 ) stop 'pb in def geometry function'
  do i = 1, lx
    x = real(i-1,dp)
    y = a*sin(2._dp*pi*x/Lx) + b + midly
    ninty = nint(y)
    mirror = midly
    do j = midly, ly
      if( j < ninty ) then
        inside(i,j,:) = fluid
      else
        inside(i,j,:) = solid
      end if
      inside(i, mirror, :) = inside(i, j, :)
      mirror = mirror - 1
    end do
  end do

END SUBROUTINE CONSTRUCT_SINUSOIDAL_WALLS_2D




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
