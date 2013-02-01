! this subroutine prints a file name 'filename' using x, y, z, array(x,y,z)

subroutine print4darray(lx,ly,lz,array,filename)

  use precision_kinds, only: dp, i2b

  implicit none
  integer(kind=i2b), intent(in) :: lx, ly, lz ! number of points in each direction
  character(*), intent(in) :: filename ! filename of .cube file. For example : "density.cube"
  real(kind=dp), intent(in), dimension ( lx, ly, lz ) :: array ! array printed in .cube file
  integer(kind=i2b) :: i, j, k ! dummy


  open(99,file=filename)

  do k=1,lz
    do j=1,ly
      do i=1,lx
        write(99,*) i,j,k,array(i,j,k)
      end do
    end do
  end do

  close(99)

end subroutine print4darray
