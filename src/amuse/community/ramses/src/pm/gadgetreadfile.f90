MODULE gadgetreadfilemod
!
! Routines to read L-Gadget2 particle data files in Fortran
! Adapted from http://astro.dur.ac.uk/~jch/password_pages/index.html
!
! CHARACTER(LEN=*) :: basefile - the base snapshot name. It may get a file id added
!
!
!
!
! Data type corresponding to gadget file header 
  TYPE gadgetheadertype
     INTEGER*4, DIMENSION(6) :: npart
     REAL*8, DIMENSION(6) :: mass
     REAL*8 :: time
     REAL*8 :: redshift
     INTEGER*4 :: flag_sfr
     INTEGER*4 :: flag_feedback
     INTEGER*4, DIMENSION(6) :: nparttotal
     INTEGER*4 :: flag_cooling
     INTEGER*4 :: numfiles
     REAL*8 :: boxsize
     REAL*8 :: omega0
     REAL*8 :: omegalambda
     REAL*8 :: hubbleparam
     INTEGER*4 :: flag_stellarage
     INTEGER*4 :: flag_metals
     INTEGER*4, DIMENSION(6)  :: totalhighword
     INTEGER*4 :: flag_entropy_instead_u
     INTEGER*4 :: flag_doubleprecision
     INTEGER*4 :: flag_ic_info
     REAL*4 :: lpt_scalingfactor
     CHARACTER, DIMENSION(48) :: unused
  END TYPE gadgetheadertype

CONTAINS

! ---------------------------------------------------------------------------

  SUBROUTINE gadgetreadheader(basename,ifile, header, ok)
!
! Read and return the gadget file header for the specified file
!
    IMPLICIT NONE
! Input parameters
    CHARACTER(LEN=*), INTENT(IN) :: basename
    INTEGER, INTENT(IN):: ifile
! Header to return
    TYPE (gadgetheadertype), INTENT(OUT) :: header
    logical, INTENT(OUT)::ok
! Internal variables
    CHARACTER(LEN=256) :: filename
    CHARACTER(LEN=4) :: fileno

    filename = TRIM(basename)
    INQUIRE(file=filename, exist=ok)
    if (.not.ok) then
       !     Generate the number to go on the end of the filename
       IF(ifile.LT.10)THEN
          WRITE(fileno,'(".",1i1.1)')ifile
       ELSE IF (ifile.LT.100)THEN
          WRITE(fileno,'(".",1i2.2)')ifile
       ELSE
          WRITE(fileno,'(".",1i3.3)')ifile
       END IF
       filename = TRIM(basename) // fileno
       INQUIRE(file=filename, exist=ok)
       if(.not.ok) then
          write(*,*) 'No file '//basename//' or '//filename
          RETURN
       end if
    end if
       
    OPEN(unit=1,file=filename,status='old',action='read',form='unformatted')
    ! Byte swapping doesn't work if you just do READ(1)header
    READ(1)header%npart,header%mass,header%time,header%redshift, &
         header%flag_sfr,header%flag_feedback,header%nparttotal, &
         header%flag_cooling,header%numfiles,header%boxsize, &
         header%omega0,header%omegalambda,header%hubbleparam, &
         header%flag_stellarage,header%flag_metals,header%totalhighword, &
         header%flag_entropy_instead_u, header%flag_doubleprecision, &
         header%flag_ic_info, header%lpt_scalingfactor
    CLOSE(1)

  END SUBROUTINE gadgetreadheader

! ---------------------------------------------------------------------------

  SUBROUTINE gadgetreadfile(basename,ifile,header,pos,vel,id)
!
! Read and return all data from the specified file. Output arrays must
! already be allocated. Use readheader to get particle numbers to do this.
!
    IMPLICIT NONE
! Input parameters
    CHARACTER(LEN=*), INTENT(IN) :: basename
    INTEGER, INTENT(IN) :: ifile
! Header and hash table to return
    TYPE (gadgetheadertype) :: header
! Particle data
    REAL, DIMENSION(3,*) :: pos,vel
    INTEGER*4, DIMENSION(*) :: id
! Internal variables
    CHARACTER(LEN=256) :: filename
    CHARACTER(LEN=4) :: fileno
    INTEGER :: np
    logical::ok

    !     Generate the number to go on the end of the filename
    IF(ifile.LT.10)THEN
       WRITE(fileno,'(".",1i1.1)')ifile
    ELSE IF (ifile.LT.100)THEN
       WRITE(fileno,'(".",1i2.2)')ifile
    ELSE
       WRITE(fileno,'(".",1i3.3)')ifile
    END IF

    filename = TRIM(basename) // fileno

    INQUIRE(file=filename, exist=ok)

    if(.not.ok) then
        write(*,*) 'No file '//filename
        RETURN
    end if
    OPEN(unit=1,file=filename,status='old',action='read',form='unformatted')
    ! Byte swapping doesn't appear to work if you just do READ(1)header
    READ(1)header%npart,header%mass,header%time,header%redshift, &
         header%flag_sfr,header%flag_feedback,header%nparttotal, &
         header%flag_cooling,header%numfiles,header%boxsize, &
         header%omega0,header%omegalambda,header%hubbleparam, &
         header%flag_stellarage,header%flag_metals,header%totalhighword, &
         header%flag_entropy_instead_u, header%flag_doubleprecision, &
         header%flag_ic_info, header%lpt_scalingfactor
    np=header%npart(2)
    READ(1)pos(1:3,1:np)
    READ(1)vel(1:3,1:np)
    READ(1)id(1:np)
    CLOSE(1)

  END SUBROUTINE gadgetreadfile

! ---------------------------------------------------------------------------

  SUBROUTINE gadgetwritefile(basename,ifile,header,pos,vel,id)
!
! Read and return all data from the specified file. Output arrays must
! already be allocated. Use readheader to get particle numbers to do this.
!
    IMPLICIT NONE
! Input parameters
    CHARACTER(LEN=*), INTENT(IN) :: basename
    INTEGER, INTENT(IN) :: ifile
! Header and hash table to return
    TYPE (gadgetheadertype) :: header
! Particle data
    REAL, DIMENSION(3,*) :: pos,vel
    INTEGER*4, DIMENSION(*) :: id
! Internal variables
    CHARACTER(LEN=256) :: filename
    CHARACTER(LEN=4) :: fileno
    INTEGER :: np
    logical::ok


    !     Generate the number to go on the end of the filename
    IF(ifile.LT.10)THEN
       WRITE(fileno,'(".",1i1.1)')ifile
    ELSE IF (ifile.LT.100)THEN
       WRITE(fileno,'(".",1i2.2)')ifile
    ELSE
       WRITE(fileno,'(".",1i3.3)')ifile
    END IF

    filename = TRIM(basename) // fileno

    OPEN(unit=1,file=filename,status='unknown',action='write',form='unformatted')
    WRITE(1)header%npart,header%mass,header%time,header%redshift, &
         header%flag_sfr,header%flag_feedback,header%nparttotal, &
         header%flag_cooling,header%numfiles,header%boxsize, &
         header%omega0,header%omegalambda,header%hubbleparam, &
         header%flag_stellarage,header%flag_metals,header%totalhighword, &
         header%flag_entropy_instead_u, header%flag_doubleprecision, &
         header%flag_ic_info, header%lpt_scalingfactor, header%unused
    np=header%npart(2)
    WRITE(1)pos(1:3,1:np)
    WRITE(1)vel(1:3,1:np)
    WRITE(1)id(1:np)

    CLOSE(1)
    END SUBROUTINE gadgetwritefile
END MODULE gadgetreadfilemod

