module velocity_profiles

use precision_kinds
use system, only: supercell, node
use constants, only : x, y, z
implicit none

contains

subroutine write_flux_field
        integer(i2b) :: i, j, k, lx, ly, lz
        real(dp) :: vx, vy, vz
        lx = supercell%geometry%dimensions%indiceMax(x)
        ly = supercell%geometry%dimensions%indiceMax(y)
        lz = supercell%geometry%dimensions%indiceMax(z)
            
            open(10,file="output/FluxField.dat")
            do i= 1,lx 
                do j= 1,ly
                    do k= 1,lz
                        vx = node(i,j,k)%solventFlux(x)!/supercell%node(i,j,k)%solventDensity
                        vy = node(i,j,k)%solventFlux(y)!/supercell%node(i,j,k)%solventDensity
                        vz = node(i,j,k)%solventFlux(z)!/supercell%node(i,j,k)%solventDensity
                        write(10,*)i, j, k, vx, vy, vz !sum(n(i,j,k,:))/(lx*ly*lz)
                    end do
                end do
            end do
            close(10)
        end subroutine

subroutine write_density_field
        integer(i2b) :: i, j, k, lx, ly, lz
        
        lx = supercell%geometry%dimensions%indiceMax(x)
        ly = supercell%geometry%dimensions%indiceMax(y)
        lz = supercell%geometry%dimensions%indiceMax(z)
        
            open(10,file="output/DensityField.dat")
            do i= 1,lx 
                do j= 1,ly
                    do k= 1,lz
                        write(10,*)i, j, k, node(i,j,k)%solventdensity
                    end do
                end do
            end do
            close(10)
 end subroutine
 
 subroutine read_flux_field
 use read_file
 implicit none
  character(len=300) :: filename
  double precision, dimension(:,:), allocatable :: tab
  integer :: tabsize(2), i
  filename = "output/FluxField.dat"
  call read_tabwnl(filename, tab, 6)
  tabsize = shape(tab)

  do i=1,tabsize(1)
     node(nint(tab(i,1)),nint(tab(i,2)),nint(tab(i,3)))%solventFlux(x) = tab(i,4)
     node(nint(tab(i,1)),nint(tab(i,2)),nint(tab(i,3)))%solventFlux(y) = tab(i,5)
     node(nint(tab(i,1)),nint(tab(i,2)),nint(tab(i,3)))%solventFlux(z) = tab(i,6)
   end do
 
 end subroutine
 
 subroutine read_density_field
 use read_file
 implicit none
  character(len=300) :: filename
  double precision, dimension(:,:), allocatable :: tab
  integer :: tabsize(2), i
  filename = "output/DensityField.dat"
  call read_tabwnl(filename, tab, 4)
  tabsize = shape(tab)
  
  do i=1,tabsize(1)
     node(nint(tab(i,1)),nint(tab(i,2)),nint(tab(i,3)))%solventDensity = tab(i,4)
  end do
 
 end subroutine


end module
