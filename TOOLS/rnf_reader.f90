PROGRAM rnf_reader
  !!======================================================================
  !!                     ***  PROGRAM  rnf_reader  ***
  !!=====================================================================
  !!  ** Purpose :  Read runoff file for a configuration
  !!
  !!  ** Method  : Use Dai and Trenberth data base. From a mask file for
  !!               river mouth (prepared with BMGTOOLS), each river mouth
  !!               area is computed, and the monthly discharge of the river
  !!               is spreaded over the respective area.
  !!               The coding rules for this program are somewhat modified in
  !!               order to follow the variables name of Dai and Trenberth
  !!
  !! History : 1.0  : 11/2014  : J.M. Molines : Original code
  !!                
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------


  !!----------------------------------------------------------------------
  !! DATATOOLS , MEOM 2014
  !! $Id$
  !!----------------------------------------------------------------------
   USE netcdf

   IMPLICIT NONE
   
   TYPE river                  ! Data structure to reflect original data base
     INTEGER      :: station
     INTEGER      :: id
     REAL(KIND=4) :: lon
     REAL(KIND=4) :: lat
     REAL(KIND=4) :: lon_mou
     REAL(KIND=4) :: lat_mou
     REAL(KIND=4) :: area_stn
     REAL(KIND=4) :: area_mou
     REAL(KIND=4) :: vol_stn
     REAL(KIND=4) :: ratio_m2s
     REAL(KIND=4) :: xnyr
     INTEGER      :: yrb
     INTEGER      :: yre
     REAL(KIND=4) :: elev
     CHARACTER(LEN=30) :: ct_name
     CHARACTER(LEN=30) :: cn_name
     CHARACTER(LEN=30) :: riv_name
     CHARACTER(LEN=30) :: stn_name
     CHARACTER(LEN=30) :: ocn_name
     REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: flow
     REAL(KIND=4), DIMENSION(12)             :: monthly_flow
   END TYPE river

   TYPE(river), DIMENSION(:), ALLOCATABLE :: runoff

   INTEGER :: js, jt, jy, jmon
   INTEGER :: ijs, ij, iye, iyb, iy, ny
   INTEGER :: ncid, ierr, idt, idstn, id
   INTEGER :: npt, npstn
   INTEGER, DIMENSION(:), ALLOCATABLE :: itime
   INTEGER, DIMENSION(12) ::  ncount_month

   ! NATL60
!   REAL(KIND=4) :: rlonmin = -82 , rlonmax = 9.375   ! define window for current model
!   REAL(KIND=4) :: rlatmin =  25 , rlatmax = 66      !  
   ! eNATL60
!  REAL(KIND=4) :: rlonmin = -100 , rlonmax = 50     ! define window for current model
!  REAL(KIND=4) :: rlatmin =    6 , rlatmax = 70     !  
   ! eNATL36X
   REAL(KIND=4) :: rlonmin = -110 , rlonmax = 50     ! define window for current model
   REAL(KIND=4) :: rlatmin =   -15 , rlatmax = 80     !  
   ! global
!  REAL(KIND=4) :: rlonmin = -180 , rlonmax = 180     ! define window for current model
!  REAL(KIND=4) :: rlatmin =    -90 , rlatmax = 90     !  
   REAL(KIND=4) :: spval                             ! fill value for flow
   REAL(KIND=4) :: zv                                ! current value

   REAL(KIND=8), DIMENSION(12) :: dsum

!  CHARACTER(LEN=80)                            :: cf_in='coastal-stns-Vol-monthly.updated-oct2007.nc'
   CHARACTER(LEN=80)                            :: cf_in='coastal-stns-Vol-monthly.updated-Aug2014.nc4'
   CHARACTER(LEN=30), DIMENSION(:), ALLOCATABLE :: cname
!------------------------------------------------------------------------
   ! Open file and read dimensions of data set
   ierr = NF90_OPEN(cf_in, NF90_NOWRITE, ncid)
   ierr = NF90_INQ_DIMID( ncid,'time'   ,idt   ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, idt  , len= npt  ) 
   ierr = NF90_INQ_DIMID( ncid,'station',idstn ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, idstn, len= npstn) 

   ALLOCATE ( runoff(npstn) )
   ALLOCATE ( itime(npt) )
   ierr = NF90_INQ_VARID (ncid, 'time', id) ; ierr = NF90_GET_VAR(ncid, id, itime )

   ! Character variables are archived as arrays of chars. Need to read them as a whole and copy
   ! into the structure
   ALLOCATE ( cname(npstn) )
   ierr = NF90_INQ_VARID (ncid, 'ct_name' , id) ; ierr = NF90_GET_VAR(ncid, id, cname ) ; runoff(:)%ct_name=cname
   ierr = NF90_INQ_VARID (ncid, 'cn_name' , id) ; ierr = NF90_GET_VAR(ncid, id, cname ) ; runoff(:)%cn_name=cname
   ierr = NF90_INQ_VARID (ncid, 'riv_name', id) ; ierr = NF90_GET_VAR(ncid, id, cname ) ; runoff(:)%riv_name=cname
   ierr = NF90_INQ_VARID (ncid, 'stn_name', id) ; ierr = NF90_GET_VAR(ncid, id, cname ) ; runoff(:)%stn_name=cname
   ierr = NF90_INQ_VARID (ncid, 'ocn_name', id) ; ierr = NF90_GET_VAR(ncid, id, cname ) ; runoff(:)%ocn_name=cname
   DEALLOCATE (cname)
 

   DO js = 1,npstn
    ALLOCATE ( runoff(js)%flow(npt) )
    ierr = NF90_INQ_VARID( ncid,'station' , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%station , start=(/js/)) !, count=(/1/)   ) 
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 1 ; endif
    ierr = NF90_INQ_VARID( ncid,'id'      , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%id      , start=(/js/)) !, count=(/1/)   ) 
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 2 ; endif
    ierr = NF90_INQ_VARID( ncid,'lon'     , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%lon     , start=(/js/)) !, count=(/1/)   ) 
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 3 ; endif
    ierr = NF90_INQ_VARID( ncid,'lat'     , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%lat     , start=(/js/)) !, count=(/1/)   ) 
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 4 ; endif
    ierr = NF90_INQ_VARID( ncid,'lon_mou' , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%lon_mou , start=(/js/)) !, count=(/1/)   ) 
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 5 ; endif
    ierr = NF90_INQ_VARID( ncid,'lat_mou' , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%lat_mou , start=(/js/)) !, count=(/1/)   ) 
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 6 ; endif
    ierr = NF90_INQ_VARID( ncid,'vol_stn' , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%vol_stn , start=(/js/)) !, count=(/1/)   ) 
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 7 ; endif
    ierr = NF90_INQ_VARID( ncid,'ratio_m2s', id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%ratio_m2s, start=(/js/)) !, count=(/1/) ) 
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 8 ; endif
    ierr = NF90_INQ_VARID( ncid,'xnyr',     id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%xnyr    , start=(/js/)) !, count=(/1/)   ) 
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 9 ; endif
    ierr = NF90_INQ_VARID( ncid,'yrb',      id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%yrb     , start=(/js/)) !, count=(/1/)   ) 
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 10 ; endif
    ierr = NF90_INQ_VARID( ncid,'yre',      id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%yre     , start=(/js/)) !, count=(/1/)   ) 
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 11 ; endif
    ierr = NF90_INQ_VARID( ncid,'elev',     id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%elev    , start=(/js/)) !, count=(/1/)   ) 
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 12 ; endif
    ierr = NF90_INQ_VARID( ncid,'FLOW',     id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%flow(:),  start = (/ js,1 /), count=(/ 1,npt /) )
     if ( ierr /= NF90_NOERR )  then ; print *, NF90_STRERROR(ierr) ; stop 18 ; endif
    IF ( js == 1 ) ierr = NF90_GET_ATT (ncid, id, '_FillValue', spval )
   ENDDO
  ierr = NF90_CLOSE(ncid)

   ! screen data ( for information)
  print *,'         No   m2s_ratio        lonm             latm             Vol(km3/yr)    CT CN OCN River_Name      Station_Name'

   ! ouput tab separated results
   DO js = 1, npstn
    IF ( (rlonmin < runoff(js)%lon_mou ) .AND. ( runoff(js)%lon_mou < rlonmax ) .AND. &
       & (rlatmin < runoff(js)%lat_mou ) .AND. ( runoff(js)%lat_mou < rlatmax ) ) THEN
      PRINT *, js,char(9), runoff(js)%ratio_m2s,char(9), runoff(js)%lon_mou,char(9), &
             runoff(js)%lat_mou,char(9), runoff(js)%vol_stn,char(9),                 &
             TRIM(runoff(js)%ct_name), char(9), TRIM(runoff(js)%cn_name),char(9),     &
             TRIM(runoff(js)%ocn_name),char(9), TRIM(runoff(js)%riv_name),char(9),    &
             TRIM(runoff(js)%stn_name)
    ENDIF
   ENDDO

! Compute monthly mean for all station
  
  DO js = 1, npstn
    runoff(js)%monthly_flow(:) = spval
    iyb=runoff(js)%yrb ; iye = runoff(js)%yre
    IF ( ( iyb == -999) .OR. (iye == -999) ) CYCLE
    ny=iye-iyb+1
    iy = (iyb- 1900)*12 +1
    ncount_month (:) = 0
    ! cumulate on dble precision
    dsum(:) = 0.d0
    DO jy = iyb, iye
     DO jmon=1,12
         ! index of point for jy, jmon
         ij= (jy-iyb)*12 + iy +jmon -1
         zv=runoff(js)%flow(ij)
         IF ( zv /= spval ) THEN
           dsum(jmon) = dsum(jmon) + DBLE(zv)
           ncount_month(jmon) = ncount_month(jmon) +1
         ENDIF
     ENDDO
    ENDDO

    ! compute  monthly mean 
    DO jmon = 1, 12
        IF (ncount_month(jmon) == 0 ) THEN 
!          print *, ' PB with station ', js, TRIM(runoff(js)%riv_name),' ', TRIM(runoff(js)%ct_name), ' ',TRIM(runoff(js)%cn_name), ' ',TRIM(runoff(js)%ocn_name),' month ', jmon, iyb, iye
        ELSE
          runoff(js)%monthly_flow(jmon) = dsum(jmon)/ncount_month(jmon)
        ENDIF
    ENDDO
!         print *,  TRIM(runoff(js)%riv_name),SUM(runoff(js)%monthly_flow(:))/12., runoff(js)%vol_stn
  ENDDO
  

END PROGRAM rnf_reader

