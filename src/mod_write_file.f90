module write_file

implicit none 

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Write output tab in a file.dat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine write_tab(filename, tab, nbline)
    implicit none
    double precision, dimension(:,:) :: tab
    integer :: tabsize(2), i, j, nbline, step
    character(len=300) :: filename
    
    j = 0
    tabsize = shape(tab)
    call execute_command_line('touch '//trim(filename))
    open(unit=2, file=filename)
    if (tabsize(1) < nbline .or. nbline == 0) then
        !write(2,*) '#', tabsize
         do i = 1,tabsize(1)
              write(2,*)tab(i,1:tabsize(2))
         end do
    elseif (tabsize(1) > nbline) then
        step = floor(real(tabsize(1))/real(nbline))
        write(2,*) '#', nbline, tabsize(2)
        do i = step, tabsize(1), step 
             write(2,*)tab(i,1:tabsize(2))
        end do
    endif
    close(2)
  end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Write Geom in vmd file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine write_vmd(filename, geom, option)
    implicit none
    integer :: VMDsolid = 4 ! pink for solid
    integer :: VMDfluid = 1 ! white for fluid
    integer, dimension(:,:,:) :: geom
    character(len=300) :: filename
    integer :: i, j, k, lx, ly, lz, temp(3), nbsolid, option
    double precision :: factor
                

    factor = 1.5d0
    temp = shape(geom)
    lx = temp(1)
    ly = temp(2)
    lz = temp(3)
    nbsolid = SUM(geom(:,:,:))

    call execute_command_line('touch '//trim(filename))
    open(452,file=filename)

    write(452,*)'PRIMVEC'
    write(452,*) REAL(lx), 0.0, 0.0
    write(452,*) 0.0, REAL(ly), 0.0
    write(452,*) 0.0, 0.0, REAL(lz)
    write(452,*)
    write(452,*)'PRIMCOORD'
    !write(452,*) lx*ly*lz, 1
    
    
    select case(option)
    case(0)
        write(452,*) lx*ly*lz, 1
        DO i=1,lx
            DO j=1,ly
            DO k=1,lz
                    write(452,*)INT((VMDsolid - VMDfluid)*geom(i,j,k) + VMDfluid), factor*REAL(i), factor*REAL(j), factor*REAL(k)
            END DO
          END DO
        END DO

    case(1)
        write(452,*) lx*ly*lz - nbsolid, 1
        DO i=1,lx
            DO j=1,ly
            DO k=1,lz
                    if (geom(i,j,k) == 0) write(452,*) 4 , factor*REAL(i), factor*REAL(j), factor*REAL(k)
            END DO
          END DO
        END DO
    
    case(2)
        write(452,*) nbsolid, 1
        DO i=1,lx
            DO j=1,ly
            DO k=1,lz
                    if (geom(i,j,k) == 1) write(452,*) 4 , factor*REAL(i), factor*REAL(j), factor*REAL(k)
            END DO
          END DO
        END DO
    end select

   CLOSE(452)
 END SUBROUTINE


end module
