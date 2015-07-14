subroutine velocity_profiles (time)
  use precision_kinds
  use system, only: supercell, node
  use constants, only : x, y, z
  implicit none
  integer(i2b), intent(in) :: time
  call VectorField

contains

  subroutine VectorField
    implicit none
    integer(i2b) :: i, j, k
    real(sp) :: vx, vy, vz
    open(10,file="output/velocityField.dat")
    do i= supercell%geometry%dimensions%indiceMin(x), supercell%geometry%dimensions%indiceMax(x)
      do j= supercell%geometry%dimensions%indiceMin(y), supercell%geometry%dimensions%indiceMax(y)
        do k= supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
          vx = REAL(node(i,j,k)%solventFlux(x),sp)!/node(i,j,k)%solventDensity
          vy = REAL(node(i,j,k)%solventFlux(y),sp)!/node(i,j,k)%solventDensity
          vz = REAL(node(i,j,k)%solventFlux(z),sp)!/node(i,j,k)%solventDensity
          write(10,*) i, j, k, vx, vy, vz
        end do
      end do
    end do
    close(10)
  end subroutine

end subroutine
