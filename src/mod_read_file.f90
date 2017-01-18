module read_file
implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Read tab of real values with unknown number of lines and ncol columns
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine read_tabwnl(filename, tab, ncol)
  implicit none
  character(len=300) :: filename
  double precision, dimension(:,:), allocatable, intent(inout) :: tab
  integer :: i, ios, l, tabsize(2), ncol
  ios = 0
  open(unit=1, file=trim(filename))
     l = 0
     do while(ios == 0)    
        l = l+1           
        read(1,*,iostat=ios)
        !print*, l, ios           
     end do
     l = l-1  !la dernière ligne est vide: on ne veut pas la lire
     !print *, 'l = ', l
     allocate(tab(1:l,1:ncol))
   close(1)

   open(unit=2, file=trim(filename))
     do i = 1,l 
           read(2,*) tab(i,1:ncol)
     end do
   close(2)

  end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Read tab of integer values with unknown number of lines and ncol columns
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine read_tabwnl_int(filename, tab, ncol)
  implicit none
  character(len=300) :: filename
  integer, dimension(:,:), allocatable, intent(inout) :: tab
  integer :: i, ios, l, tabsize(2), ncol
  ios = 0
  open(unit=1, file=trim(filename))
     l = 0
     do while(ios == 0)    
        l = l+1           
        read(1,*,iostat=ios)
        !print*, l, ios           
     end do
     l = l-1  !la dernière ligne est vide: on ne veut pas la lire
     !print *, 'l = ', l
     allocate(tab(1:l,1:ncol))
   close(1)

   open(unit=2, file=trim(filename))
     do i = 1,l 
           read(2,*) tab(i,1:ncol)
     end do
   close(2)

  end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Read last line of ncol tab in filename
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine read_tab_lastline(filename, vect, ncol)
  implicit none
  character(len=300) :: filename
  double precision, dimension(:), allocatable, intent(inout) :: vect
  integer :: ios, ncol
  ios = 0
  allocate(vect(1:ncol))

  open(unit=1, file=trim(filename))
     do while(ios == 0)    
        read(1,*,iostat=ios) vect(1:ncol)       
     end do   
   close(1)

  end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! read table of double in file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! the size of the table must be written at the first line of the file

  subroutine read_tab(filename,tab)
    implicit none
    character(len=300) :: filename
    character(len=1) :: flag
    double precision, dimension(:,:), allocatable, intent(inout) :: tab
    double precision, dimension(:,:), allocatable :: tabtemp
    integer :: nline, ncol, i, ios
    
    ios = 4
    open(unit=1, file=filename, iostat=ios)
    read(1,*) flag, nline, ncol
    allocate(tabtemp(1:nline,1:ncol))
    
    loop: do i = 1,nline
        read(1,*, iostat=ios) tabtemp(i,1:ncol)
        if (ios > 0) then 
            EXIT loop
      end if
    end do loop
    close(1)

    if (ios > 0) then
        allocate(tab(1:i-1,1:ncol))
        tab(1:i-1,:) = tabtemp(1:i-1,:)
    else
        allocate(tab(1:nline,1:ncol))
        tab(:,:) = tabtemp(:,:)
    end if
    
  end subroutine 
end module
