MODULE GEOMETRY
  use precision_kinds, only: dp, i2b
  use system, only: lx, ly, lz, fluid, solid
  implicit none
  private
  public construct_wall, construct_cylinder, construct_bcc, construct_fcc

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CONSTRUCT_WALL(inside)
  implicit none
  integer(i2b), dimension(:,:,:), intent(out) :: inside
  inside = fluid
  inside( :, :, lbound(inside,3)) = solid ! the lower bound of the thrid dimension of inside is solid
  inside( :, :, ubound(inside,3)) = solid ! so is the upper bound
END SUBROUTINE CONSTRUCT_WALL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CONSTRUCT_CYLINDER(inside)
  implicit none
  integer(i2b), dimension(:,:,:), intent(out) :: inside
  real(dp) :: radius ! radius of cylinder
  real(dp), dimension(2) :: rnode, rorigin ! coordinates of each node and center of cylinder in x,y coordinates
  integer(i2b) :: i, j ! dummy

  if( lx /= ly) stop 'wall=2 is for cylinders, which should have same lx and ly'
  rorigin = (/ real(lx,dp)/2.0_dp, real(ly,dp)/2.0_dp /)
  radius = real(lx,dp)/2.0_dp

  do concurrent (i=1:lx, j=1:ly)
    rnode = (/ i, j /)
    if( norm2(rnode-rorigin) >= radius ) then
      inside(i,j,:) = solid
    else
      inside(i,j,:) = fluid
    end if
  end do
END SUBROUTINE CONSTRUCT_CYLINDER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CONSTRUCT_BCC(inside)
  implicit none
  integer(i2b), dimension(:,:,:), intent(out) :: inside
  integer(i2b) :: i,j,k ! dummy
  if( lx /= ly .or. lx /= lz ) stop 'with wall = 3, i.e. cfc cell, the supercell should be cubic with lx=ly=lz'
  inside = fluid
  do concurrent( i=1:lx, j=1:ly, k=1:lz )
    if( is_in_bcc_solid_sphere(i,j,k) ) then
      inside(i,j,k) = solid
    else
      inside(i,j,k) = fluid
    end if
  end do
  contains
    PURE LOGICAL FUNCTION IS_IN_BCC_SOLID_SPHERE(i,j,k)
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
        is_in_bcc_solid_sphere = .true.
      else
        is_in_bcc_solid_sphere = .false.
      end if
    END FUNCTION IS_IN_BCC_SOLID_SPHERE
END SUBROUTINE CONSTRUCT_BCC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CONSTRUCT_FCC(inside)
  implicit none
  integer(i2b), dimension(:,:,:), intent(out) :: inside
  integer(i2b) :: i,j,k ! dummy
  if( lx /= ly .or. lx /= lz ) stop 'with wall = 3, i.e. cfc cell, the supercell should be cubic with lx=ly=lz'
  inside = fluid
  do concurrent( i=1:lx, j=1:ly, k=1:lz )
    if( is_in_fcc_solid_sphere(i,j,k) ) then
      inside(i,j,k) = solid
    else
      inside(i,j,k) = fluid
    end if
  end do
  contains
    PURE LOGICAL FUNCTION IS_IN_FCC_SOLID_SPHERE(i,j,k)
      implicit none
      integer(i2b), intent(in) :: i,j,k
      real(dp), dimension(14) :: distances
      distances(1) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/0.0_dp, 0.0_dp, 0.0_dp/) )
      distances(2) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp), 0.0_dp, 0.0_dp/) )
      distances(3) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/0.0_dp, real(ly,dp), 0.0_dp/) )
      distances(4) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/0.0_dp, 0.0_dp, real(lz,dp)/) )
      distances(5) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp), real(ly,dp), 0.0_dp/) )
      distances(6) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp), 0.0_dp, real(lz,dp)/) )
      distances(7) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/0.0_dp, real(ly,dp), real(lz,dp)/) )
      distances(8) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp), real(ly,dp), real(lz,dp)/) )
      distances(9) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp)/2.0_dp, real(ly,dp)/2.0_dp, 0.0_dp/) )
      distances(10) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp)/2.0_dp, 0.0_dp, real(lz,dp)/2.0_dp/) )
      distances(11) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/0.0_dp, real(ly,dp)/2.0_dp, real(lz,dp)/2.0_dp/) )
      distances(12) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp)/2.0_dp, real(ly,dp)/2.0_dp, real(lz,dp)/) )
      distances(13) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp)/2.0_dp, real(ly,dp), real(lz,dp)/2.0_dp/) )
      distances(14) = norm2( (/real(i,dp),real(j,dp),real(k,dp)/) - (/real(lx,dp), real(ly,dp)/2.0_dp, real(lz,dp)/2.0_dp/) )
      ! put distances to all of corners of the cube and to its center (=9 points) in an array
      if ( any(distances <= real(lx,dp)*sqrt(2.0_dp)/4.0_dp) ) then
        is_in_fcc_solid_sphere = .true.
      else
        is_in_fcc_solid_sphere = .false.
      end if
    END FUNCTION IS_IN_FCC_SOLID_SPHERE
END SUBROUTINE CONSTRUCT_FCC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE GEOMETRY
