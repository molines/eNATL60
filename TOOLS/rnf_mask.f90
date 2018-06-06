PROGRAM rnf_mask
  !!======================================================================
  !!                     ***  PROGRAM  rnf_mask  ***
  !!=====================================================================
  !!  ** Purpose : mask the input file according to the bathymetry given
  !!               as arguments.
  !!
  !!  ** Method  : Read both file and set data file to zero where bathy is 0
  !!
  !! History :  1.0  : 06/2018  : J.M. Molines : 
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
  INTEGER :: npiglo, npjglo, nx, ny, nxx, nyy
  INTEGER :: ncid, id, ierr, ncidd, iddat
  INTEGER :: narg, ijarg

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: bathy  ! on eDomain
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdata  ! on eDomain

  CHARACTER(LEN=80) :: cf_bathy
  CHARACTER(LEN=80) :: cf_data
  CHARACTER(LEN=80) :: cf_out='none'
  CHARACTER(LEN=80) :: cv_data  = 'Bathymetry'
  CHARACTER(LEN=80) :: cv_bathy = 'Bathymetry'
  CHARACTER(LEN=80) :: cldum
  !!------------------------------------------------------------------------

  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  rnf_mask -d DATA-file -b BATHY-file [-vd VAR-data] [-vb VAR-bathy]'
     PRINT *,'          [-o FILE-out] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Create masked copy of the input data file'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -d DATA-file  : give the name of the data file'
     PRINT *,'       -b BATHY-file : give the name of the bathymetry '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -vd VAR-data  : give the name of the variable in data file'
     PRINT *,'            default is ',TRIM(cv_data)
     PRINT *,'       -vb VAR-bathy  : give the name of the variable in bathymetry file'
     PRINT *,'            default is ',TRIM(cv_bathy)
     PRINT *,'       -o FILE-out: specify name of output file, instead of <BATHY-new>.wrk'
     PRINT *,'      '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : <DATA-file>.wrk or  <FILE-OUT> '
     PRINT *,'         variables : ',TRIM(cv_data),' or VAR-data if specified'
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
        ! options
     CASE ( '-vd'  ) ; CALL getarg(ijarg, cv_data  ) ; ijarg=ijarg+1
     CASE ( '-vb'  ) ; CALL getarg(ijarg, cv_bathy ) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out   ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  IF ( TRIM(cf_out) == 'none' ) THEN
     cf_out=TRIM(cf_data)//".wrk"
  ENDIF

  CALL system ( "cp "//TRIM(cf_data)//" "//TRIM(cf_out))
  ! work on cf_out directly
  cf_data=cf_out
  PRINT *,' WORKING ON ',TRIM(cf_data)

  ! Read size of edomain :
  PRINT *,' OPENING '//TRIM(cf_data)
  ierr= NF90_OPEN(cf_data, NF90_WRITE, ncid)
  ierr= NF90_INQ_DIMID(ncid,'x',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid, id, len=npiglo )
  ierr= NF90_INQ_DIMID(ncid,'y',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid, id, len=npjglo )

  PRINT *,'     NPIGLO = ', npiglo
  PRINT *,'     NPJGLO = ', npjglo

  ! allocate bathy and mask array
  ALLOCATE ( bathy(npiglo, npjglo), rdata(npiglo, npjglo) )
  ierr = NF90_INQ_VARID(ncid,cv_data,iddat)
  ierr = NF90_GET_VAR(ncid,iddat,rdata)

  ! Open and read data file
  ierr = NF90_OPEN(cf_bathy,NF90_NOWRITE,ncidd)
  ierr = NF90_INQ_DIMID(ncidd,'x',id) ; ierr=NF90_INQUIRE_DIMENSION(ncidd, id, len=nxx )
  ierr = NF90_INQ_DIMID(ncidd,'y',id) ; ierr=NF90_INQUIRE_DIMENSION(ncidd, id, len=nyy )

  IF ( nxx /= npiglo .OR. nyy /= npjglo ) THEN
     PRINT *, ' ERROR : mismatch in specified window size, incoherent with data file.'
     STOP
  ENDIF
  ierr = NF90_INQ_VARID(ncidd,cv_bathy,id)
  ierr = NF90_GET_VAR(ncidd,id,bathy)

  WHERE (bathy == 0. ) rdata = 0.

  ! write back modified data to file
  ierr=NF90_PUT_VAR(ncid,iddat,rdata)
  ierr=NF90_CLOSE(ncid)

END PROGRAM rnf_mask

