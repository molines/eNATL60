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
  INTEGER :: npiglo, npjglo, nx, ny, nrim, nxx, nyy
  INTEGER :: ncid, id, ierr, ncidd, idbat
  INTEGER :: narg, ijarg
  INTEGER :: iimin, iimax, ijmin, ijmax
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: mask        ! on eDomain, inner part on Domain

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: bathy  ! on eDomain
  REAL(KIND=4)                              :: rdepmin=3 !m

  CHARACTER(LEN=80) :: cf_bathy='eNATL60_BATHY_GEBCO_2014_2D_raw_v2.nc'
  CHARACTER(LEN=80) :: cf_data
  CHARACTER(LEN=80) :: cf_out='none'
  CHARACTER(LEN=80) :: cv_data
  CHARACTER(LEN=80) :: cldum
  !!------------------------------------------------------------------------
  !! initialize default corner values
  iimin=1053 ; iimax=6474
  ijmin=1276 ; ijmax=4729
  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  data_extend_zone -b BATHY-new -d DATA-old -v VAR-old ...'
     PRINT *,'          [-o FILE-out] [-w imin imax jmin jmax ]'
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
     PRINT *,'       -o FILE-out: specify name of output file, instead of <BATHY-new>.wrk'
     PRINT *,'      '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : Bathymetry !! '
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
        ;            ; CALL getarg(ijarg, cldum  )   ; ijarg=ijarg+1 ; READ(cldum,*) ijmax
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out )   
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

  IF ( TRIM(cf_out) == 'none' ) THEN
     cf_out=TRIM(cf_bathy)//".wrk"
  ENDIF

  CALL system ( "cp "//TRIM(cf_bathy)//" "//TRIM(cf_out))
  ! work on cf_out directly
  cf_bathy=cf_out
  PRINT *,' WORKING ON ',TRIM(cf_bathy)

  ! Read size of edomain :
  PRINT *,' OPENING '//TRIM(cf_bathy)
  ierr= NF90_OPEN(cf_bathy, NF90_WRITE, ncid)
  ierr= NF90_INQ_DIMID(ncid,'x',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid, id, len=npiglo )
  ierr= NF90_INQ_DIMID(ncid,'y',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid, id, len=npjglo )

  PRINT *,'     NPIGLO = ', npiglo
  PRINT *,'     NPJGLO = ', npjglo
  PRINT *,'         NX = ', nx
  PRINT *,'         NY = ', ny

  ! allocate bathy and mask array
  ALLOCATE ( bathy(npiglo, npjglo), mask(npiglo, npjglo) )
  ierr = NF90_INQ_VARID(ncid,'Bathymetry',idbat)
  ierr = NF90_GET_VAR(ncid,idbat,bathy)

  ! read Domain mask
  ! build the  mask
  mask=1
  WHERE ( bathy /= 0. ) bathy=1.
  ! Open and read data file
  ierr = NF90_OPEN(cf_data,NF90_NOWRITE,ncidd)
  ierr = NF90_INQ_DIMID(ncidd,'x',id) ; ierr=NF90_INQUIRE_DIMENSION(ncidd, id, len=nxx )
  ierr = NF90_INQ_DIMID(ncidd,'y',id) ; ierr=NF90_INQUIRE_DIMENSION(ncidd, id, len=nyy )

  IF ( nxx /= nx .OR. nyy /= ny ) THEN
     PRINT *, ' ERROR : mismatch in specified window size, incoherent with data file.'
     STOP
  ENDIF
  ierr = NF90_INQ_VARID(ncidd,cv_data,id)
  ierr = NF90_GET_VAR(ncidd,id,bathy(iimin:iimax,ijmin:ijmax) )
 

  ! write back modifies bathy to file
  ierr=NF90_PUT_VAR(ncid,idbat,bathy)
  ierr=NF90_CLOSE(ncid)

END PROGRAM data_extend_zone

