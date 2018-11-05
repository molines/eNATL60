PROGRAM rnf_merge
  !!======================================================================
  !!                     ***  PROGRAM  rnf_merge  ***
  !!=====================================================================
  !!  ** Purpose : Replace runoff in an existing file by new runoff from
  !!               a new file
  !!
  !!  ** Method  : Read and replace
  !!
  !! History :  1.0  : 11/2018  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------


  !!----------------------------------------------------------------------
  !! RUNOFF , MEOM 2018
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER(KIND=4) :: ierr, ncid_old, ncid_new, id
  INTEGER(KIND=4) :: narg, ijarg
  INTEGER(KIND=4) :: npiglo, npjglo, npt

  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: runoff_old, runoff_new 

  CHARACTER(LEN=80) :: cf_old
  CHARACTER(LEN=80) :: cf_new
  CHARACTER(LEN=80) :: cf_out
  CHARACTER(LEN=80) :: cv_old='sorunoff'
  CHARACTER(LEN=80) :: cv_new='sorunoff'
  CHARACTER(LEN=80) :: cldum

  !========================================================================
  narg=iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  rnf_merge -old OLD-runoff -new NEW-runoff  -o OUT-runoff '
     PRINT *,'          [ -vo OLD-var ] [-vn NEW-var] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'         Create OUT-file as a runoff file where NEW-runoff replace'
     PRINT *,'        old-runoff where they are present.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -old OLD-runoff : name of the old runoff file '
     PRINT *,'       -new NEW-runoff : name of the new runoff file '
     PRINT *,'       -o OUT-runoff : name of the merged runoff file '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        [-vo OLD-var ] : give name of OLD runoff var [',TRIM(cv_old),']'
     PRINT *,'        [-vn NEW-var ] : give name of NEW runoff var [',TRIM(cv_new),']'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cv_new),' (    )'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      ' 
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-old'   ) ; CALL getarg(ijarg, cf_old ) ; ijarg=ijarg+1
     CASE ( '-new'   ) ; CALL getarg(ijarg, cf_new ) ; ijarg=ijarg+1
     CASE ( '-o'     ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1
        ! option
     CASE ( '-vo'    ) ; CALL getarg(ijarg, cv_old ) ; ijarg=ijarg+1
     CASE ( '-vn'    ) ; CALL getarg(ijarg, cv_new ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO
  IF ( cf_out == cf_new ) THEN
      PRINT *,' NEW and OUT should be different !'
      STOP
  ENDIF

  ! copy cf_new to cf_out and work on cf_out 
  CALL system ('cp '//TRIM(cf_new)//' '//TRIM(cf_out) )

  ! Read dimension in cv_old
  ierr = NF90_OPEN(cf_old, NF90_NOWRITE, ncid_old )
  ierr = NF90_INQ_DIMID(ncid_old,'x',            id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid_old,id,len=npiglo)
  ierr = NF90_INQ_DIMID(ncid_old,'y',            id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid_old,id,len=npjglo)
  ierr = NF90_INQ_DIMID(ncid_old,'time_counter', id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid_old,id,len=npt)
  PRINT *,'  NPIGLO = ', npiglo
  PRINT *,'  NPJGLO = ', npjglo

  ! allocate arrays
  ALLOCATE ( runoff_old(npiglo,npjglo,npt) , runoff_new(npiglo,npjglo,npt) )

  ierr = NF90_INQ_VARID(ncid_old,cv_old,id ) ; ierr = NF90_GET_VAR(ncid_old, id, runoff_old )
  ierr = NF90_CLOSE(ncid_old)

  ! Now read new file
  ierr = NF90_OPEN(cf_out, NF90_WRITE, ncid_new )
  ierr = NF90_INQ_VARID(ncid_new,cv_new,id ) 
  ierr = NF90_GET_VAR(ncid_new, id, runoff_new )
  ! now merge
  WHERE ( runoff_new == 0. ) runoff_new=runoff_old
  ierr = NF90_PUT_VAR(ncid_new, id, runoff_new )

  ! deal now with socoeff variable 
  ierr = NF90_INQ_VARID(ncid_new,'socoefr',id ) 
  runoff_old=0.
  WHERE ( runoff_new(:,:,1) /= 0. ) runoff_old(:,:,1)=0.5
  ierr = NF90_PUT_VAR(ncid_new, id, runoff_old(:,:,1) )

  
  ierr = NF90_CLOSE(ncid_new) 
  
END PROGRAM rnf_merge
