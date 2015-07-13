module io
    use precision_kinds
    use system, only: supercell, node
    use constants, only: x, y, z
    use input , only : input_line
    implicit none
    contains

        subroutine manage
            call put_input_in_character_array
            call print_input_in_output_folder ! Print input parameters found to output folder
        end subroutine

        ! this subroutine prints a file name 'filename' using x, y, z, array(x,y,z)
        subroutine print4darray(lx,ly,lz,array,filename )
            integer(i2b), intent(in) :: lx, ly, lz ! number of points in each direction
            character(*), intent(in) :: filename ! filename of .cube file. For example : "density.cube"
            real(dp), intent(in), dimension ( lx, ly, lz ) :: array ! array printed in .cube file
            integer(i2b) :: i, j, k ! dummy
            open(99,file=filename)
            do k=1,lz
                do j=1,ly
                    do i=1,lx
                        write(99,*) i,j,k,array(i,j,k)
                    end do
                end do
            end do
            close(99)
        end subroutine

        subroutine print_everything_related_to_charge_equil
            use system, only: phi
            call print4darray(ubound(phi,1),ubound(phi,2),ubound(phi,3),phi,'output/phi_of_x_y_z.dat') ! print internal potential
        end subroutine

        ! In this subroutine, one prints the first lines of codes to stdout.
        subroutine print_header
            character(8)  :: date
            character(10) :: time
            call date_and_time ( DATE=date,TIME=time)
            print*,
            print*,
            print*,date(1:4),'/',date(5:6),'/',date(7:8),', ',time(1:2),'h',time(3:4),'m',time(5:6)
            print*,'=================================================='
            print*,'Laboetie, fluid dynamics for chemical applications'
            print*,'=================================================='
            print*,
        end subroutine

        ! this subroutine prints all the input parameters in output/input.out
        ! so that all files needed to understand the outputs are in the output folder.
        ! written by Maximilien Levesque, 2011, @ Ecole Normale Superieure
        subroutine print_input_in_output_folder
            use input, only: input_line
            integer(i2b) :: i ! dummy for loop
            open( unit = 99 , file='output/lb.in' ) ! open file to write in
            ! print each line of input_line()
            do i = 1 , size( input_line )
                write(99,*) input_line (i)
            end do
        end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Print an XSF file of the supercell for visualisation in VMD for instance.
        ! Type vmd --xsf output/supercell.xsf to visualise it.
        SUBROUTINE print_supercell_xsf

            USE system ,ONLY: fluid, solid, supercell
            ! Lx size of box in angstroms
            ! nb_solute_sites nombre de sites (pareil que le nombre de lignes dans format xyz)
            ! xmol is an array which contains les x de tous les sites, il a donc la taille de x_mol(nb_solute_sites)
            ! tu peux remplacer integer(kind=i2b) par integer tout court.
            INTEGER :: i, j, k, lx, ly, lz
            INTEGER, PARAMETER :: xsfUnit=5

            lx = supercell%geometry%dimensions%indiceMax(x)
            ly = supercell%geometry%dimensions%indiceMax(y)
            lz = supercell%geometry%dimensions%indiceMax(z)

            open(xsfUnit,file='output/supercell.xsf')

            100 format (xA)
            101 format (3(xF10.5))
            102 format (xI5,xI1)
            103 format (xI3,3(xxF10.5))

            write(xsfUnit,100)'# this is the specification file of the supercell'
            write(xsfUnit,100)'# lines beginning with # are commented. There cannot be comment lines within the sections'
            write(xsfUnit,100)'# XSF format specifications can be found at http://www.xcrysden.org/doc/XSF.html'
            write(xsfUnit,100)'# I strongly recommends to read this documentation.'
            write(xsfUnit,*)
            write(xsfUnit,100)'# for periodic structures one has to begin with word CRYSTAL'
            write(xsfUnit,100)'CRYSTAL'
            write(xsfUnit,100)
            write(xsfUnit,100)'# Then one needs to specify the lattice vectors'
            write(xsfUnit,100)'# specification of PRIMVEC (in ANGSTROMS) like:'
            write(xsfUnit,100)'#         ax, ay, az    (first lattice vector)'
            write(xsfUnit,100)'#         bx, by, bz    (second lattice vector)'
            write(xsfUnit,100)'#         cx, cy, cz    (third lattice vector)'
            write(xsfUnit,100)'# pay attention to vectors as they are written in horizontal way which is quite unusual'
            write(xsfUnit,100)'# for now only orthorhombic structures allowed (a/=b/=c, all angles are 90 degrees)'
            write(xsfUnit,100)'PRIMVEC'
            write(xsfUnit,101) real( [lx,0,0] ,sp)
            write(xsfUnit,101) real( [0,ly,0] ,sp)
            write(xsfUnit,101) real( [0,0,lz] ,sp)
            write(xsfUnit,*)
            write(xsfUnit,100)'# Then one needs to specify the atoms belonging to the unit cell. '
            write(xsfUnit,100)'# First number stands for number of atoms in the primitive cell (2 in this case).'
            write(xsfUnit,100)'# The second number is always 1 for PRIMCOORD coordinates.'
            write(xsfUnit,100)'# in angstroms and cartesian coordinates'
            write(xsfUnit,100)'PRIMCOORD'
            write(xsfUnit,*) lx*ly*lz, 1

            BLOCK
                REAL(sp) :: coo(3)
                INTEGER :: nature
                LOGICAL :: isInterfacial
                INTEGER, PARAMETER :: VMDpink=4, VMDgreen=45, VMDwhite=1

                DO i=1,lx
                    DO j=1,ly
                        DO k=1,lz

                            nature        = node(i,j,k)%nature
                            isInterfacial = node(i,j,k)%isInterfacial
                            coo = REAL([i-1,j-1,k-1],sp)

                            IF      ( nature == solid )                            THEN
                                WRITE(5,*) VMDpink,coo
                                ! PRINT*,coo,"solid"
                            ELSE IF ( nature == fluid .AND. isInterfacial )        THEN
                                WRITE(5,*) VMDgreen,coo
                                ! PRINT*,coo,"fluid interfacial"
                            ELSE IF ( nature == fluid .AND. (.NOT.isInterfacial) ) THEN
                                WRITE(5,*) VMDwhite,coo
                                ! PRINT*,coo,"fluid NOT interfacial"
                            ELSE
                                STOP 'While writing supercell.xsf, I found a node that has undocumented nature'
                            END IF

                        END DO
                    END DO
                END DO

            END BLOCK

            CLOSE(5)
        END SUBROUTINE

!subroutine print_vacf
!  use system, only: dp, i2b, vacf, tmax, tmom
!  use constants, only: x, y, z
!  implicit none
!  integer(kind=i2b) :: t
!  open(unit=99,file='output/vacf.dat')
!
!  if( ubound(vacf,2) > tmax-tmom+1) stop 'ubound vacfx too low !?'
!  do t= 0, tmax-tmom
!    write(99,*) t, vacf(x,t), vacf(y,t), vacf(z,t)
!  end do
!  print*,'wrote output/vacf.dat'
!  close(99)
!end subroutine print_vacf

        ! This subroutine writes an intent(in) array called 'array' to 'filename' in the file format .cube that handles 3-dimensional data.
        ! The standard .cube file can be read by vmd for instance. It is meant to vizualize periodic supercells.
        ! see for instance http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html
        ! or http://paulbourke.net/dataformats/cube/
        ! It is a very easy format to deal with because it's very natural to write and read.
        ! Inputs needed : array to write, filename, number of points in each direction, length of each supercell vector, atomic number of the
        ! atoms and their coordinates. The total number of atoms is determined as the number of different coordinates ie as size(x_mol)
        ! output : a file called 'filename' containing 3dimensional data in .cube format
        ! written by Maximilien Levesque, while in post doc at Ecole Normale Superieure, Paris in Daniel Borgis's theoretical chemistry group
        ! 20110919  Maximilien Levesque, clean version for Virginie M.

        subroutine write_to_cube_file ( Lx, Ly, Lz, array, filename )
            integer(kind=i2b) :: i, j, k!, n ! dummy variables
            integer(kind=i2b), intent(in) :: lx, ly, lz
            character(*), intent(in) :: filename ! filename of .cube file. For example : "density.cube"
            real(kind=dp), intent(in), dimension ( lx, ly, lz ) :: array ! array printed in .cube file
            real(kind=dp), parameter :: angtobohr = 1.889725989_dp ! 1 Bohr = 1.889725989 Ang. Necessary because of VMD understanding of lengths
            ! define formats for writing to file
            ! 104 format ( xI3 , xA , 3(xxF10.5) )
            102 format ( 1(xF10.5) )
            ! first open the file you want to print array in. It's a formatted file.
            open(10, file = filename , form = 'formatted' )
            write(10,*) ' CPMD CUBE FILE.' ! default text
            write(10,*) ' OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z' ! default text, for remembering
            ! write the total number of sites and a default text nobody knows it meaning
            write( 10 , * ) 1 ,' 0.0 0.0 0.0 ' ! 0 0 0 ou Lx/2 Ly/2 Lz/2 ?  Size(x_mol) is the total number of atoms
            ! write primary vectors
            write( 10 , * ) lx, lx*angtobohr, ' 0.0 0.0'
            write( 10 , * ) ly, ' 0.0 ' , ly* angtobohr , ' 0.0'
            write( 10 , * ) lz, ' 0.0 0.0 ' , lz* angtobohr
            ! write the atoms and their coordinates in Bohr
            write(10,*) '1 0.0 0.0 0.0 0.0'
            ! write .cube file. One value per line. As said earlier, run over x then y then z. The one varying the most rapidly is z.
            do i = 1 , lx
                do j = 1 , ly
                    do k = 1 , lz
                        write ( 10 , 102 ) array ( i , j , k )
                    end do
                end do
            end do
            close(10) ! close the .cube file called filename
            print*, filename,' written' ! warn user
        end subroutine

        subroutine inquireNecessaryFilesExistence
            logical :: file_exists
            inquire(file="./output/.", exist=file_exists)
            if( .not. file_exists) then
                call system('mkdir -p output') ! -p do not print error if exist and create parent directory if needed
            end if
            inquire(file="./lb.in", exist=file_exists)
            if( .not. file_exists) stop "./lb.in, i.e. the input file, does not exist. STOP."
        end subroutine

        subroutine put_input_in_character_array
            integer(i2b) :: i, j, k, n ! dummy
            integer(i2b) :: ios ! input output status of readen file
            character(len=100) :: text ! temporary input line
            character(len=100), allocatable, dimension(:) :: arraytemp  ! Temporary array to stock data for resizing input_line
            character(len=len('lb.in')), parameter :: inputFilename = "lb.in"
            integer(i2b) :: totalnumberofinputlines
            ! open the file in which you have all inputs
            open ( unit = 11 , file = inputFilename )
            ! count the total number of lines (n) in input.in
            ! a blank line counts for 0 if it is only a return to next line (\n)
            ! a line with only blanks counts for 1
            totalnumberofinputlines = 0
            do while (.true.)
                read ( 11 , '(a)' , iostat = ios ) text
                if ( ios > 0 ) then
                    write (*,*) 'error in input. critical'
                    stop
                else if ( ios < 0 ) then
                    exit ! exit do loop and stop reading file
                else
                    if ( text /= '\n' ) totalnumberofinputlines = totalnumberofinputlines + 1
                end if
            end do
            ! close file in order to open it properly from start
            close ( unit = 11 )
            ! allocate input_line so that it has the right number of lines
            allocate ( input_line (totalnumberofinputlines) )
            input_line = ' ' ! init
            ! open again the input file
            open ( unit = 11 , file = inputfilename )
            !read it and put each line in input_line
            do i = 1 , totalnumberofinputlines
                read ( 11 , '(a)'  ) text
                input_line (i) = trim ( adjustl (text) ) ! trim () removes trailing blanks while adjustl removes left blanks and put white blanks at the end
            end do
            ! clean up comments in the lines (expl: option = 3 # blabla)
            do i = 1 , totalnumberofinputlines
                do j = 1 , len(text)
                    if ( input_line (i) (j:j) == '#' ) then
                        forall ( k = j : len(text) )
                            input_line (i) (k:k) = ' '
                        end forall
                        exit
                    end if
                end do
                input_line (i) = trim ( adjustl (input_line (i) ) )
            end do
            !Delete blank lines and count the size of the smallest array containing initial data
            n=0 !init
            do i = 1 , totalnumberofinputlines
                if ( input_line (i) (1:1) /= ' ' )  then
                    input_line (n+1) = input_line(i)
                    n = n + 1
                endif
            end do
            !Resize input_line to the smallest size by using a temporary array
            arraytemp = input_line
            deallocate ( input_line )
            allocate ( input_line ( n ) )
            input_line = arraytemp ( 1 : n  )
        end subroutine put_input_in_character_array

end module
