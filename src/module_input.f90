MODULE input

  USE precision_kinds ,ONLY: i2b,dp
  IMPLICIT NONE
  CHARACTER(len=100), ALLOCATABLE, DIMENSION(:) :: input_line ! array containing all input lines
  LOGICAL :: verbose
  PRIVATE
  PUBLIC :: verbose, input_line, input_dp, input_int, input_log, input_char, n_linesInFile, deltaAbscissa, &
    input_dp2, input_dp3, input_int2, input_int3

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE FUNCTION input_dp (tag, defaultValue)
    IMPLICIT NONE
    REAL(DP) :: input_dp
    CHARACTER(*), INTENT(IN) :: tag
    REAL(DP), optional, intent(in) :: defaultValue
    INTEGER(i2b) :: i, j
    logical :: ifoundtag
    ifoundtag = .false.
    j=LEN(tag)
    DO i = 1, SIZE( input_line)
      IF( input_line( i)( 1:j) == tag  .AND. input_line(i)(j+1:j+1)==' ' ) then
        READ(input_line(i)(j+4:j+50),*) input_dp
        ifoundtag = .true.
        exit
      end if
    END DO
    if (ifoundtag .eqv. .false. .and. present(defaultValue)) input_dp = defaultValue
  END FUNCTION input_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE FUNCTION input_dp2 (tag, defaultValue)
    IMPLICIT NONE
    REAL(DP) :: input_dp2(2)
    CHARACTER(*), INTENT(IN) :: tag
    REAL(DP), optional, intent(in) :: defaultValue(2)
    INTEGER(i2b) :: i, j
    logical :: ifoundtag
    ifoundtag = .false.
    j=LEN(tag)
    DO i = 1, SIZE( input_line)
      IF( input_line( i)( 1:j) == tag  .AND. input_line(i)(j+1:j+1)==' ' ) then
        READ(input_line(i)(j+4:j+50),*) input_dp2
        ifoundtag = .true.
        exit
      end if
    END DO
    if (ifoundtag .eqv. .false. .and. present(defaultValue)) input_dp2 = defaultValue
  END FUNCTION input_dp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE FUNCTION input_dp3 (tag, defaultValue)
    IMPLICIT NONE
    REAL(DP) :: input_dp3(3)
    CHARACTER(*), INTENT(IN) :: tag
    REAL(DP), optional, intent(in) :: defaultValue(3)
    INTEGER(i2b) :: i, j
    logical :: ifoundtag
    ifoundtag = .false.
    j=LEN(tag)
    DO i = 1, SIZE( input_line)
      IF( input_line( i)( 1:j) == tag  .AND. input_line(i)(j+1:j+1)==' ' ) then
        READ(input_line(i)(j+4:j+50),*) input_dp3
        ifoundtag = .true.
        exit
      end if
    END DO
    if (ifoundtag .eqv. .false. .and. present(defaultValue)) input_dp3 = defaultValue
  END FUNCTION input_dp3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE FUNCTION input_int (tag, defaultValue)
    IMPLICIT NONE
    INTEGER(I2B) :: input_int
    CHARACTER(*), INTENT(IN) :: tag
    integer(i2b), optional, intent(in) :: defaultValue
    INTEGER(i2b) :: i, j
    logical :: ifoundtag
    ifoundtag = .false.
    j=LEN(tag)
    DO i = 1, SIZE( input_line)
      IF( input_line( i)( 1:j) == tag  .AND. input_line(i)(j+1:j+1)==' ' ) then
        READ(input_line(i)(j+4:j+50),*) input_int
        ifoundtag = .true.
        exit
      end if
    END DO
    if (ifoundtag.eqv..false. .and. present(defaultValue)) input_int = defaultValue
  END FUNCTION input_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE FUNCTION input_int2 (tag, defaultValue)
    IMPLICIT NONE
    INTEGER(I2B) :: input_int2(2)
    CHARACTER(*), INTENT(IN) :: tag
    integer(i2b), optional, intent(in) :: defaultValue(2)
    INTEGER(i2b) :: i, j
    logical :: ifoundtag
    ifoundtag = .false.
    j=LEN(tag)
    DO i = 1, SIZE( input_line)
      IF( input_line( i)( 1:j) == tag  .AND. input_line(i)(j+1:j+1)==' ' ) then
        READ(input_line(i)(j+4:j+50),*) input_int2
        ifoundtag = .true.
        exit
      end if
    END DO
    if (ifoundtag.eqv..false. .and. present(defaultValue)) input_int2 = defaultValue
  END FUNCTION input_int2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE FUNCTION input_int3 (tag, defaultValue)
    IMPLICIT NONE
    INTEGER(I2B) :: input_int3(3)
    CHARACTER(*), INTENT(IN) :: tag
    integer(i2b), optional, intent(in) :: defaultValue(3)
    INTEGER(i2b) :: i, j
    logical :: ifoundtag
    ifoundtag = .false.
    j=LEN(tag)
    DO i = 1, SIZE( input_line)
      IF( input_line( i)( 1:j) == tag  .AND. input_line(i)(j+1:j+1)==' ' ) then
        READ(input_line(i)(j+4:j+50),*) input_int3
        ifoundtag = .true.
        exit
      end if
    END DO
    if (ifoundtag.eqv..false. .and. present(defaultValue)) input_int3 = defaultValue
  END FUNCTION input_int3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION input_log (tag, defaultValue)
    IMPLICIT NONE
    logical :: input_log
    CHARACTER(*), INTENT(IN) :: tag
    logical, intent(in), optional :: defaultValue
    CHARACTER :: text
    INTEGER(i2b) :: i, j, lentag
    logical :: found
    found = .false.
    IF (tag=='point_charge_electrostatic') THEN
      STOP 'The tag point_charge_electrostatic in dft.in must be renamed direct_sum since July 27th, 2014'
    END IF
    lentag = LEN(tag)
    DO i =1, SIZE( input_line)
      IF( input_line(i)(1:lentag)==tag .AND. input_line(i)(lentag+1:lentag+1)==' ' ) then
        READ( input_line(i)(lentag+4:lentag+50) ,*) text
        found = .true.
        exit
      end if
    END DO
    if( .not.found .and. present(defaultValue) ) then
      input_log = defaultValue
      return
    else if( .not.found ) then
      print*, "I could not find keyword '", tag,"' in ./input/dft.in"
      print*, "It should have been there associated to logical T or F"
      stop
    end if
    j = 999 ! means error in reading
    IF( text(1:1) == 'T' ) j = 1 ! means true, 2 means false
    IF( text(1:1) == 't' ) j = 1
    IF( text(1:1) == 'F' ) j = 2
    IF( text(1:1) == 'f' ) j = 2
    IF( j == 999 ) THEN
      PRINT*, 'Problem with logical tag "', tag,' in ./input/dft.in'
      print*, "It is here but is not logical. I read from it: '",text
      STOP
    END IF
    IF( j == 1 ) input_log = .TRUE.
    IF( j == 2 ) input_log = .FALSE.
  END FUNCTION input_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION input_char (tag, defaultvalue)
    IMPLICIT NONE
    CHARACTER(50) :: input_char
    CHARACTER(*), INTENT(IN) :: tag
    INTEGER(i2b) :: i,j,imax,iostatint
    CHARACTER(*), INTENT(IN), OPTIONAL :: defaultValue
    j=LEN(tag)
    i=0
    imax=SIZE(input_line)
    DO i=1,imax+1
      IF (i==imax+1) THEN
        IF( PRESENT(defaultvalue) ) THEN
            input_char = defaultvalue
            RETURN
        ELSE
            PRINT*,"I didnt find keyword '",tag,"'in input file."
            ERROR STOP
        END IF
      END IF
      IF (input_line(i)(1:j)==tag .AND. input_line(i)(j+1:j+1)==' ') THEN
        READ(input_line(i)(j+4:j+50),*,IOSTAT=iostatint) input_char
        IF (iostatint/=0) THEN
          PRINT*,"I have a problem in reading input line:"
          PRINT*,TRIM(ADJUSTL(input_line(i)))
          IF (TRIM(ADJUSTL(input_line(i)(j+4:j+50)))=='') PRINT*,"I found nothing after sign ="
          STOP
        END IF
        EXIT
      END IF
    END DO
    IF (input_char(1:1)==' ') THEN
      PRINT*,"First character of ",tag," is a whitespace"
      STOP
    END IF
    IF (LEN(TRIM(ADJUSTL(input_char)))==0) THEN
      PRINT*,"Tag after ",tag," is only whitespaces."
      STOP
    END IF
  END FUNCTION input_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION n_linesInFile (filename)
    IMPLICIT NONE
    INTEGER :: n_linesInFile
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER :: ios
    OPEN (77, FILE=filename)
    n_linesInFile = 0
    DO WHILE (.true.)
      READ (77,*,IOSTAT=ios)
      IF (ios>0) THEN
        WRITE(*,*)'Error in file:',filename
        STOP
        ELSE IF (ios<0) THEN ! end of file reached
          EXIT
        ELSE
          n_linesInFile = n_linesInFile +1
        END IF
      END DO
      CLOSE (77)
  END FUNCTION n_linesInFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION deltaAbscissa (filename)
    REAL(dp) :: abscissa, previousAbscissa, ordonates, deltaAbscissa
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER(i2b) :: i, ios, n_lines
    OPEN (10, FILE=filename, IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE(*,*)"Cant open file ",filename," in FUNCTION deltaAbscissa"
    END IF

    DO i= 1, 2
      READ(10,*,IOSTAT=ios) previousAbscissa, ordonates
      IF (ios/=0) then
        PRINT*, 'Something went wrong while reading ', TRIM(ADJUSTL(filename)), ' in module_input>deltaAbscissa'
        STOP
      END IF
      READ(10,*, IOSTAT=ios) abscissa, ordonates
      IF (ios/=0) then
        PRINT*, 'Something went wrong while reading ', TRIM(ADJUSTL(filename)), ' in module_input>deltaAbscissa'
        STOP
      END IF
    END DO
    deltaAbscissa = abscissa-previousAbscissa
    CLOSE(10)
    n_lines = n_linesInFile(filename)
    OPEN (10, FILE=filename, IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE(*,*)"Cant open file ",filename," in FUNCTION deltaAbscissa"
    END IF

    READ(10,*,IOSTAT=ios) abscissa, ordonates
    IF (ios/=0) then
      PRINT*, 'Something went wrong while reading ', TRIM(ADJUSTL(filename)), ' in module_input>deltaAbscissa'
      STOP
    END IF
    DO i=1, n_lines-1
      previousAbscissa = abscissa
      READ(10,*,IOSTAT=ios) abscissa, ordonates
      IF (ios>0) then
        PRINT*, 'Something went wrong while reading ', TRIM(ADJUSTL(filename)), ' in module_input>deltaAbscissa'
        STOP
      ELSE IF (ios<0) THEN
        EXIT
      ELSE
        IF ((abscissa-previousAbscissa-deltaAbscissa)/deltaAbscissa > 1E-5) THEN
          PRINT*, abscissa, previousAbscissa,abscissa-previousAbscissa ,deltaAbscissa
          PRINT*, 'STOP. Non uniform absissa in ', filename
          STOP
        END IF
      END IF
    END DO
    CLOSE(10)
  END FUNCTION deltaAbscissa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE input