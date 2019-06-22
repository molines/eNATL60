PROGRAM  bat_mark_pool
  !!======================================================================
  !!                     ***  PROGRAM  bat_mark_pool  ***
  !!=====================================================================
  !!  ** Purpose : detection of unconnected points in a bathymetric file
  !!
  !!  ** Method  : use fill_pool2D algo
  !!
  !! History :  1.0  : 06/2019  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   fill_pool2D   : taken from CDFTOOLS ( P. Mathiot)
  !!----------------------------------------------------------------------
  USE netcdf

  IMPLICIT NONE

  INTEGER(KIND=4) :: ncid, ierr, id, npiglo, npjglo
  INTEGER(KIND=4) :: ncidg
  INTEGER(KIND=4) :: iiseed, ijseed, ji,jj
  INTEGER(KIND=4) :: narg, ijarg, iargc

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: bathy
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: bathyg

  CHARACTER(LEN=80) :: cf_bat
  CHARACTER(LEN=80) :: cf_bat_pool
  CHARACTER(LEN=80) :: cf_bat_ground
  CHARACTER(LEN=80) :: cdum

  LOGICAL :: ln_ground=.false.

  !!----------------------------------------------------------------------
  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *, 'USAGE : bat_mark_pool -b BATHY-file -s ISEED JSEED [-ground]'
     PRINT *
     PRINT *,'   PURPOSE: '
     PRINT *,'       work on a copy of BATHY_file and fill all the point connected'
     PRINT *,'       to the seed point (ISEED,JSEED) with -1. Unconnected point thus'
     PRINT *,'       are the only positive points in the resulting file. This makes'
     PRINT *,'       easier the elimination of these unconnected points.'
     PRINT *,'   '
     PRINT *,'   ARGUMENT:'
     PRINT *,'      -b BATHY-file : give the name of the bathymetric file. Variable name'
     PRINT *,'            is assumed to be Bathymetry.'
     PRINT *,'      -s ISEED JSEED : gives the I,J position of the seed point. This point'
     PRINT *,'            should be in the ocean.'
     PRINT *,'   '
     PRINT *,'   OPTIONS:'
     PRINT *,'       -ground : This option will ground unconnected ocean points (setting '
     PRINT *,'                 them to 0. Results will be output in a specific file.' 
     PRINT *,' '
     PRINT *,'   OUTPUT:'
     PRINT *,'       Information is given on the standard output.'
     PRINT *,'       <BATHY-file>.pool  contains the Bathymetry field after filling process'
     PRINT *,'       Positive values corresponds to unconnected points.... to be fixed.'
     PRINT *,'       Note that the program let the user free to fix the unconnected points.'
     PRINT *,'       The grounding of these points (setting them to 0) is not the only '
     PRINT *,'       solution ! (Adding some ocean points may allow the connection for'
     PRINT *,'       isolated pools. '
     PRINT *,'       If you decide to just ground the unconnected points,'
     PRINT *,'       use the -ground option. The resulting bathymetry will be output in the'
     PRINT *,'       <BATHY-file>.grounded.'
     STOP
  ENDIF

  ijarg=1
  DO WHILE (ijarg < narg )
     CALL getarg (ijarg, cdum) ; ijarg = ijarg + 1
     SELECT CASE ( cdum)
     CASE ( '-b' ) ; CALL getarg(ijarg, cf_bat) ; ijarg = ijarg + 1
     CASE ( '-s' ) ; CALL getarg(ijarg, cdum  ) ; ijarg = ijarg + 1  ; READ(cdum,*) iiseed
        ; CALL getarg(ijarg, cdum  ) ; ijarg = ijarg + 1  ; READ(cdum,*) ijseed
     CASE ( '-ground' ) ; ln_ground=.true.
     CASE DEFAULT 
        PRINT *, 'Unkown option : ',TRIM(cdum) ; STOP
     END SELECT
  ENDDO

  cf_bat_pool=TRIM(cf_bat)//'.pool'
  CALL SYSTEM ('cp '//TRIM(cf_bat)//' '//TRIM(cf_bat_pool))
  IF (ln_ground )  THEN
     ALLOCATE (bathyg(npiglo,npjglo) )
     cf_bat_ground=TRIM(cf_bat)//'.grounded'
     CALL SYSTEM ('cp '//TRIM(cf_bat)//' '//TRIM(cf_bat_ground))
  ENDIF

     ierr=NF90_OPEN(cf_bat_pool,NF90_WRITE,ncid)
     IF (ln_ground ) ierr=NF90_OPEN(cf_bat_ground,NF90_WRITE,ncidg)
     ierr=NF90_INQ_DIMID(ncid,'x',id) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npiglo)
     ierr=NF90_INQ_DIMID(ncid,'y',id) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npjglo)

     ALLOCATE (bathy(npiglo,npjglo) )
     ierr = NF90_INQ_VARID(ncid,'Bathymetry',id)
     ierr = NF90_GET_VAR(ncid, id, bathy)
     PRINT *,NF90_STRERROR(ierr)
     IF (ln_ground ) bathyg(:,:) = bathy(:,:)

     CALL FillPool2D(iiseed, ijseed, bathy, -1. )
     print *, count( (bathy > 0 ) )
     DO ji=1, npiglo
        DO jj=1,npjglo
           IF ( bathy(ji,jj) > 0. ) THEN
              PRINT *, ji,jj,bathy(ji,jj)
           ENDIF
        ENDDO
     ENDDO
     ierr = NF90_PUT_VAR(ncid, id, bathy)
     ierr = NF90_CLOSE(ncid)
     IF (ln_ground) THEN
        WHERE (bathy > 0 ) bathyg=0.
        ierr =NF90_PUT_VAR(ncidg,id,bathyg)
        ierr = NF90_CLOSE(ncidg)
     ENDIF


   CONTAINS

     SUBROUTINE FillPool2D(kiseed, kjseed, pdta, pifill)
       !!---------------------------------------------------------------------
       !!                  ***  ROUTINE FillPool2D  ***
       !!  
       !! ** Purpose :  Replace all area surrounding by mask value by kifill value
       !!  
       !! ** Method  :  flood fill algorithm
       !!  
       !!----------------------------------------------------------------------
       INTEGER(KIND=4),              INTENT(in)    :: kiseed, kjseed
       REAL(KIND=4),                 INTENT(in)    :: pifill         ! pool value
       REAL(KIND=4), DIMENSION(:,:), INTENT(inout) :: pdta           ! mask

       INTEGER :: ik                       ! number of point change
       INTEGER :: ip                       ! size of the pile
       INTEGER :: ji, jj                   ! loop index

       INTEGER :: iip1, iim1, ii, ij       ! working integer
       INTEGER :: ijp1, ijm1               ! working integer
       INTEGER :: ipiglo, ipjglo           ! size of the domain, infered from kdta size

       INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: ipile    ! pile variable
       REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdata    ! new data
       !!----------------------------------------------------------------------
       ! infer domain size from input array
       ipiglo = SIZE(pdta,1)
       ipjglo = SIZE(pdta,2)

       ! allocate variable
       ALLOCATE(ipile(2*ipiglo*ipjglo,2))
       ALLOCATE(rdata(ipiglo,ipjglo))

       ! initialise variables
       rdata=pdta
       ipile(:,:)=0
       ipile(1,:)=[kiseed,kjseed]
       ip=1; ik=0

       ! loop until the pile size is 0 or if the pool is larger than the critical size
       DO WHILE ( ip /= 0 ) ! .AND. ik < 600000);
          ik=ik+1
          ii=ipile(ip,1); ij=ipile(ip,2)
          IF ( MOD(ik, 10000) == 0 ) PRINT *, 'IP =', ip, ik, ii,ij

          ! update bathy and update pile size
          rdata(ii,ij) =pifill
          ipile(ip,:)  =[0,0]; ip=ip-1

          ! check neighbour cells and update pile ( assume E-W periodicity )
          iip1=ii+1; IF ( iip1 == ipiglo+1 ) iip1=2
          iim1=ii-1; IF ( iim1 == 0        ) iim1=ipiglo-1

          ijp1=ij+1; IF ( ijp1 == ipjglo+1 ) ijp1=ipjglo
          ijm1=ij-1; IF ( ijm1 == 0        ) ijm1=1

          IF (rdata(ii, ijp1) > 0 ) THEN
             ip=ip+1; ipile(ip,:)=[ii  ,ijp1]
          END IF
          IF (rdata(ii, ijm1) > 0 ) THEN
             ip=ip+1; ipile(ip,:)=[ii  ,ijm1]
          END IF
          IF (rdata(iip1, ij) > 0 ) THEN
             ip=ip+1; ipile(ip,:)=[iip1,ij  ]
          END IF
          IF (rdata(iim1, ij) > 0 ) THEN
             ip=ip+1; ipile(ip,:)=[iim1,ij  ]
          END IF

       END DO
       pdta=rdata;

       DEALLOCATE(ipile); DEALLOCATE(rdata)

     END SUBROUTINE FillPool2D





   END PROGRAM bat_mark_pool
