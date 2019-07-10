PROGRAM rnf_12_36
  !!======================================================================
  !!                     ***  PROGRAM  rnf_12_36  ***
  !!=====================================================================
  !!  ** Purpose : transform rnf file from 1/12 to 1/36
  !!
  !!  ** Method  : Each pixel of 1/12 is copied on 9 pixel of 1/36
  !!               
  !!
  !! History :  1.0  : 06/2019  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------

  USE netcdf
  !!----------------------------------------------------------------------
  !!  eNATL60 Tools .
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: ji,jj, jt
  INTEGER :: ii36, ij36
  INTEGER :: npiglo12, npjglo12, nx12, ny12, npt12
  INTEGER :: npiglo36, npjglo36, nx36, ny36
  INTEGER :: ncid12, id, ierr
  INTEGER :: ncid36, idx, idy, idt, id36, idcof
  INTEGER :: narg, ijarg

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdata12  ! on eDomain12
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdata36  ! on eDomain36

  CHARACTER(LEN=80) :: cf_in12
  CHARACTER(LEN=80) :: cf_out36
  CHARACTER(LEN=80) :: cldum
  !!------------------------------------------------------------------------

  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  rnf_12_36 -i FILE.12  -o FILE.36'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Create a refined copy of FILE.12'
     PRINT *,'        Each pixel of FILE.12 is splitted into 9 sub-pixel having'
     PRINT *,'        The same value.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -i FILE.12 : name for 1/12 file'
     PRINT *,'       -o FILE.36 : name for 1/36 output file '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'      '
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-i'   ) ; CALL getarg(ijarg, cf_in12   ) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out36  ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ! Read size of edomain :
  PRINT *,' OPENING '//TRIM(cf_in12)
  ierr= NF90_OPEN(cf_in12, NF90_WRITE, ncid12)
  ierr= NF90_INQ_DIMID(ncid12,'x',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid12, id, len=npiglo12 )
  ierr= NF90_INQ_DIMID(ncid12,'y',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid12, id, len=npjglo12 )
  ierr= NF90_INQ_DIMID(ncid12,'time_counter',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid12, id, len=npt12 )

  PRINT *,'     NPIGLO = ', npiglo12
  PRINT *,'     NPJGLO = ', npjglo12
  PRINT *,'     TIME_COUNTER = ', npt12

  npiglo36=npiglo12*3
  npjglo36=npjglo12*3

! Create output file
  ierr = NF90_CREATE(cf_out36,NF90_NETCDF4,ncid36)
  
  ierr = NF90_DEF_DIM(ncid36,'x',npiglo36,idx)
  ierr = NF90_DEF_DIM(ncid36,'y',npjglo36,idy)
  ierr = NF90_DEF_DIM(ncid36,'time_counter',NF90_UNLIMITED,idt)

  ierr = NF90_DEF_VAR(ncid36,'runoff',NF90_FLOAT,(/idx,idy,idt/), id36)
  ierr = NF90_DEF_VAR(ncid36,'socoefr',NF90_FLOAT,(/idx,idy/), idcof)

  ! allocate bathy and mask array
  ALLOCATE ( rdata12(npiglo12, npjglo12))
  ALLOCATE ( rdata36(npiglo36, npjglo36))
  ierr = NF90_INQ_VARID(ncid12,'runoff',id)
  DO jt = 1, npt12
  ierr = NF90_GET_VAR(ncid12,id,rdata12,start=(/1,1,jt/), count=(/npiglo12,npjglo12,1/))
  DO jj=1,npjglo12
    DO ji=1,npiglo12
       ii36=(ji-1)*3+2
       ij36=(jj-1)*3+2
       rdata36(ii36-1:ii36+1,ij36-1:ij36+1) = rdata12(ji,jj)
    ENDDO
  ENDDO
  ierr = NF90_PUT_VAR(ncid36,id36,rdata36,start=(/1,1,jt/), count=(/npiglo36,npjglo36,1/) )
  ENDDO
  ierr = NF90_INQ_VARID(ncid12,'socoefr',id)
  ierr = NF90_GET_VAR(ncid12,id,rdata12,start=(/1,1/), count=(/npiglo12,npjglo12/))

  DO jj=1,npjglo12
    DO ji=1,npiglo12
       ii36=(ji-1)*3+2
       ij36=(jj-1)*3+2
       rdata36(ii36-1:ii36+1,ij36-1:ij36+1) = rdata12(ji,jj)
    ENDDO
  ENDDO
  ierr = NF90_PUT_VAR(ncid36,idcof,rdata36,start=(/1,1/), count=(/npiglo36,npjglo36/) )
  ierr = NF90_CLOSE(ncid12)
  ierr = NF90_CLOSE(ncid36)
  

END PROGRAM rnf_12_36

