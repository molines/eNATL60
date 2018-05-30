PROGRAM data_extend_zone
  !!======================================================================
  !!                     ***  PROGRAM  data_extend_zone  ***
  !!=====================================================================
  !!  ** Purpose : Create an extended file and copy the original values 
  !!               on the corresponding domain.
  !!
  !!  ** Method  : Use the eXtended bathymetry to set the extended mask. Then
  !!               replace copy data on the respective place from old file.
  !!               Values outside the old domain are set to some default values.
  !!
  !! History :  1.0  : 05/2018  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------

  USE netcdf
  !!----------------------------------------------------------------------
  !!  eNATL60 Tools .
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: ji,jj
  INTEGER :: npiglo, npjglo, nx, ny, nrim
  INTEGER :: ncid, id, ierr, ncidm
  INTEGER :: narg, ijarg
  INTEGER :: iimin, iimax, ijmin, ijmax
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: mask        ! on eDomain, inner part on Domain

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: bathy  ! on eDomain
  REAL(KIND=4)                              :: rdepmin=3 !m

  CHARACTER(LEN=80) :: cf_bathy='eNATL60_BATHY_GEBCO_2014_2D_raw_v2.nc'
  CHARACTER(LEN=80) :: cf_data
  CHARACTER(LEN=80) :: cv_data
  CHARACTER(LEN=80) :: cldum
  !!------------------------------------------------------------------------
  !! initialize default corner values
  iimin=1053 ; iimax=6474
  ijmin=1276 ; ijmax=4729
  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  data_extend_zone -b BATHY-new -d DATA-old -v VAR-old ...'
     PRINT *,'          [-w imin imax jmin jmax ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Create a file of the size of the new bathymetry, with old values copied'
     PRINT *,'         from old data file' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -b BATHY-new : give the name of the new bathymetry '
     PRINT *,'       -d DATA-old  : give the name of the old data file'
     PRINT *,'       -v VAR-old  : give the name of the variable to copy in the old data file'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -w imin imax jmin jmax : specify the position of old grid with respect ' 
     PRINT *,'          to the new one. Default is for NATL60 into eNATL60'
     PRINT *,'      '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cv_out),' (    )'
     PRINT *,'      '
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-b'   ) ; CALL getarg(ijarg, cf_bathy ) ; ijarg=ijarg+1
     CASE ( '-d'   ) ; CALL getarg(ijarg, cf_data  ) ; ijarg=ijarg+1
     CASE ( '-v'   ) ; CALL getarg(ijarg, cv_data  ) ; ijarg=ijarg+1
     ! options
     CASE ( '-w'   ) ; CALL getarg(ijarg, cldum  )   ; ijarg=ijarg+1 ; READ(cldum,*) iimin
        ;            ; CALL getarg(ijarg, cldum  )   ; ijarg=ijarg+1 ; READ(cldum,*) iimax
        ;            ; CALL getarg(ijarg, cldum  )   ; ijarg=ijarg+1 ; READ(cldum,*) ijmin
        ;            ; CALL getarg(ijarg, cldum  )   ; ijarg=ijarg+1 ; READ(cldum,*) ijmmax
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ! initialise corner value : ie index of Domain corners in eDomain reference

  !      4 --------  3
  !      |           |
  !      |           |
  !      1 --------  2
  !  icorner=(/1053,6474,6474,1053/)   ! v0-v1
  !  jcorner=(/ 641, 641,4094,4094/)   ! v0-v1

  nx=iimax - iimin +1
  ny=ijmax - ijmin +1

  CALL system ( "cp "//TRIM(cf_bathy)//" "//TRIM(cf_bathy)//".wrk" )
  cf_bathy=TRIM(cf_bathy)//".wrk"

  ! Read size of edomain :
  PRINT *,' OPENING '//TRIM(cf_bathy)
  ierr=NF90_OPEN(cf_bathy, NF90_WRITE, ncid)
  ierr= NF90_INQ_DIMID(ncid,'x',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid, id, len=npiglo )
  ierr= NF90_INQ_DIMID(ncid,'y',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid, id, len=npjglo )

  PRINT *,'     NPIGLO = ', npiglo
  PRINT *,'     NPJGLO = ', npjglo
  PRINT *,'         NX = ', nx
  PRINT *,'         NY = ', ny

  ! allocate bathy and mask array
  ALLOCATE ( bathy(npiglo, npjglo), mask(npiglo, npjglo) )

  ! read Domain mask
  mask=2
  ierr=NF90_OPEN(cf_mask,NF90_NOWRITE, ncidm)
  ierr=NF90_INQ_VARID(ncidm,'tmask',id)
  ierr=NF90_GET_VAR(ncidm, id, mask(iimin:iimax,ijmin:ijmax), start=(/1,1,1,1/),count=(/nx,ny,1,1/) )
  ierr=NF90_CLOSE(ncidm)

  ! modify mask : put 2 where we do not want to change the raw bathy 
  !  (1) on the rim ( 3*15 +1 points SOUTH, EAST and North --no WEST -- )
  nrim=3*15+1
  mask(icorner(1):icorner(2)        , jcorner(1):jcorner(1)+nrim   ) = 2 
  mask(icorner(2):icorner(2)-nrim:-1, jcorner(2):jcorner(3)        ) = 2
  mask(icorner(4):icorner(2)        , jcorner(4):jcorner(4)-nrim:-1) = 2 
  ! (2) Do not mask hudson bay
  mask(1:1855,jcorner(4)-1432:jcorner(4)) = 2 

  ! read eBathy :
  ierr = NF90_INQ_VARID(ncid,'Bathymetry',id)
  ierr = NF90_GET_VAR(ncid,id,bathy)

  ! mask bathy
  WHERE ( mask == 0 ) bathy=0.
  ! put minimum depth ( 3m) where mask is 1 and bathy=0
  WHERE ( bathy == 0. .AND. mask == 1 ) bathy=rdepmin

  ! write back modifies bathy to file
  ierr=NF90_PUT_VAR(ncid,id,bathy)
  ierr=NF90_CLOSE(ncid)

END PROGRAM data_extend_zone

