PROGRAM rnf_compute_area
  !!======================================================================
  !!                     ***  PROGRAM  rnf_compute_area  ***
  !!=====================================================================
  !!  ** Purpose : Compute monthly runoff file from Dai Trenberth data base
  !!
  !!  ** Method  : Read river_mouth file, model metrics and Dai Trenberth.
  !!           
  !!
  !! History : 1.0  : 11/2014  : J.M. Molines : Original code
  !!                  06/2018  : J.M. Molines : adaption for eNatl60
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------


  !!----------------------------------------------------------------------
  !! DATA_TOOLS , MEOM 2014
  !! $Id$
  !! Copyright (c) 2014, J.-M. Molines
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER :: js, jm
  INTEGER :: npi, npj, np
  INTEGER :: ncid, id, ierr
  INTEGER :: npstn
  INTEGER :: idx, idy, idt, idlon, idlat, idtime, idcoef, idrnf
  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: irivmsk, itmsk, itmp

  REAL(KIND=4),    DIMENSION(:,:),   ALLOCATABLE :: rlon, rlat
  REAL(KIND=4),    DIMENSION(:,:),   ALLOCATABLE :: socoefr
  REAL(KIND=4),    DIMENSION(:,:,:), ALLOCATABLE ::  sorunoff

  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: e1t, e2t
  REAL(KIND=8)      :: darea, dconvcoef, dconv2mm

  CHARACTER(LEN=80) :: cf_coo ='coordinates.nc'
  CHARACTER(LEN=80) :: cf_riv ='river_mouth.nc'
  CHARACTER(LEN=80) :: cf_in='coastal-stns-Vol-monthly.updated-oct2007.nc'
  CHARACTER(LEN=80) :: cf_out='runoff.nc'
  CHARACTER(LEN=80) :: cv_e1t = 'e1t'
  CHARACTER(LEN=80) :: cv_e2t = 'e2t'
  CHARACTER(LEN=80) :: cv_gla = 'glamt'
  CHARACTER(LEN=80) :: cv_gph = 'gphit'
  CHARACTER(LEN=80) :: cv_rivmsk = 'Bathymetry'   !  Stupid BMGTOOLS requirement !

  ! Dai and trenberth structure
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
     CHARACTER(LEN=30) :: ocn_name
     REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: flow
     REAL(KIND=4), DIMENSION(12)             :: monthly_flow
  END TYPE river

  TYPE(river), DIMENSION(:), ALLOCATABLE :: runoff
  !!----------------------------------------------------------------------
  ! read e1t,  e2t in coordinates.nc
  ierr = NF90_OPEN(cf_coo, NF90_NOWRITE, ncid)
  ierr = NF90_INQ_DIMID(ncid, 'x', id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npi)
  ierr = NF90_INQ_DIMID(ncid, 'y', id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npj)
  PRINT  *, ' Size of Coordinate file :',npi,npj

  ALLOCATE( e1t(npi,npj), e2t(npi,npj), rlon(npi,npj), rlat(npi,npj) )
  ierr = NF90_INQ_VARID( ncid, cv_e1t, id ) ; ierr = NF90_GET_VAR(ncid, id, e1t,  start=(/1,1,1/) )
  ierr = NF90_INQ_VARID( ncid, cv_e2t, id ) ; ierr = NF90_GET_VAR(ncid, id, e2t,  start=(/1,1,1/) )
  ierr = NF90_INQ_VARID( ncid, cv_gla, id ) ; ierr = NF90_GET_VAR(ncid, id, rlon, start=(/1,1,1/) )
  ierr = NF90_INQ_VARID( ncid, cv_gph, id ) ; ierr = NF90_GET_VAR(ncid, id, rlat, start=(/1,1,1/) )
  PRINT * ,e1t( 500,500), e2t( 500,500)
  PRINT * ,e1t( 1000,1000), e2t( 1000,1000)

  ierr = NF90_CLOSE(ncid) 

  ! read river mask
  ierr = NF90_OPEN(cf_riv, NF90_NOWRITE, ncid)
  ierr = NF90_INQ_DIMID(ncid, 'x', id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npi)
  ierr = NF90_INQ_DIMID(ncid, 'y', id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npj)

  PRINT  *, ' Size of River mask file  :',npi,npj

  ! read river data base
  CALL GetDaiTrenberth ()
  PRINT *, ' NPSTN', npstn

  ! JMM : check size !
  ALLOCATE( irivmsk(npi,npj), itmsk(npi,npj), itmp(npi,npj) )
  ALLOCATE (socoefr(npi,npj), sorunoff(npi,npj,12) )
  socoefr(:,:) = 0.
  sorunoff(:,:,:) = 0.

  ierr = NF90_INQ_VARID( ncid, cv_rivmsk, id ) ; ierr = NF90_GET_VAR(ncid, id, irivmsk, start=(/1,1,1/) )
  ierr = NF90_CLOSE(ncid)

  dconvcoef = 10.**9/365./86400.d0   ! transform km3/yr to m3/s
  dconv2mm  = 1000.d0 * 86400.
  itmsk=0
  WHERE ( irivmsk > 1 ) itmsk = irivmsk
  WHERE ( irivmsk > 1 ) socoefr=0.5
  DO js = 1, npstn
     np = COUNT ( (itmsk == js ) )
     IF ( np /= 0 ) THEN
        itmp(:,:)=0
        WHERE ( itmsk == js ) itmp = 1
        darea = SUM( e1t *e2t, (itmp == 1) )
        ! compute runoff ( taking into account some particular case ( St Lawrence, Ottawa, Saguenay for instance)
        SELECT CASE ( js )
        CASE ( 16 ) !  ( St Lawrence, Ottawa, Saguenay Ottawa St Maurice Outardes Manicuagan )
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + &
                &    runoff( 59)%vol_stn * runoff( 59)%ratio_m2s + &   !  Ottawa
                &    runoff(130)%vol_stn * runoff(130)%ratio_m2s + &   !  St Maurice
                &    runoff( 66)%vol_stn * runoff( 66)%ratio_m2s + &   !  Saguenay
                &    runoff(230)%vol_stn * runoff(230)%ratio_m2s + &   !  Outardes
                &    runoff(101)%vol_stn * runoff(101)%ratio_m2s       !  Manicuagan
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + &
                &            runoff( 59)%monthly_flow(:) * runoff( 59)%ratio_m2s + &   !  Ottawa
                &            runoff(130)%monthly_flow(:) * runoff(130)%ratio_m2s + &   !  St Maurice
                &            runoff( 66)%monthly_flow(:) * runoff( 66)%ratio_m2s + &   !  Saguenay
                &            runoff(230)%monthly_flow(:) * runoff(230)%ratio_m2s + &   !  Outardes
                &            runoff(101)%monthly_flow(:) * runoff(101)%ratio_m2s       !  Manicuagan
           runoff( js)%ratio_m2s = 1.  
        CASE ( 29 ) ! Esequibo Cuyuni
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Esequibo
                &    runoff( 93)%vol_stn * runoff( 93)%ratio_m2s                !  Cuyuni
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + &
                &    runoff( 93)%monthly_flow(:) * runoff( 93)%ratio_m2s        !  Cuyuni
           runoff( js)%ratio_m2s = 1.
        CASE ( 43 ) !  Usumasinta, Grivalva, Rapido de Sama
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Usumasinta (MX)
                &    runoff(133)%vol_stn * runoff(133)%ratio_m2s            + & ! Rapido de Sama
                &    runoff(201)%vol_stn * runoff(201)%ratio_m2s                ! Grivalva
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Usumasinta (MX)
                &    runoff(133)%monthly_flow(:) * runoff(133)%ratio_m2s    +&  ! Rapido de Sama
                &    runoff(201)%monthly_flow(:) * runoff(201)%ratio_m2s        ! Grivalva
           runoff( js)%ratio_m2s = 1.
        CASE ( 67)  !  Alabama + Tombigbee
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Alabama
                &    runoff(123)%vol_stn * runoff(123)%ratio_m2s                ! Tombigbee
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & !Alabama
                &    runoff(123)%monthly_flow(:) * runoff(123)%ratio_m2s        ! Tombigbee
           runoff( js)%ratio_m2s = 1.
        CASE (288)  ! Sabines + Noches
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Sabines
                &    runoff(317)%vol_stn * runoff(317)%ratio_m2s                ! Noches
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Sabines
                &    runoff(317)%monthly_flow(:) * runoff(317)%ratio_m2s        ! Noches
           runoff( js)%ratio_m2s = 1.
        CASE (432 ) ! Ouergha + Sehou
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Ouergah
                &    runoff(498)%vol_stn * runoff(498)%ratio_m2s                ! Sehou
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Ouergha
                &    runoff(498)%monthly_flow(:) * runoff(498)%ratio_m2s        ! Sehou
           runoff( js)%ratio_m2s = 1.
        END SELECT
        PRINT *, js, np , darea/1.e6, 'km^2', runoff(js)%vol_stn * runoff(js)%ratio_m2s/darea * dconvcoef*dconv2mm, TRIM(runoff(js)%riv_name), SUM(runoff(js)%monthly_flow(:))/12.*runoff(js)%ratio_m2s/darea*1000.*86400. ! mm/day

        DO jm=1,12
           sorunoff(:,:,jm) = sorunoff(:,:,jm) + itmp(:,:)*runoff(js)%monthly_flow(jm)*runoff(js)%ratio_m2s/darea*1000.
        ENDDO
     ENDIF
  ENDDO

  ! create runoff output file
  ierr = NF90_CREATE(cf_out,or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid)
  ierr = NF90_DEF_DIM(ncid, 'x', npi, idx )
  ierr = NF90_DEF_DIM(ncid, 'y', npj, idy )
  ierr = NF90_DEF_DIM(ncid, 'time_counter', NF90_UNLIMITED, idt )

  ierr = NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, (/ idx, idy /), idlon )
  ierr = NF90_PUT_ATT(ncid, idlon,'units','degrees_east')
  ierr = NF90_PUT_ATT(ncid, idlon,'valid_min',MINVAL(rlon) )
  ierr = NF90_PUT_ATT(ncid, idlon,'valid_max',MAXVAL(rlon) )
  ierr = NF90_PUT_ATT(ncid, idlon,'long_name','Longitude')
  ierr = NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, (/ idx, idy /), idlat )
  ierr = NF90_PUT_ATT(ncid, idlat,'units','degrees_east')
  ierr = NF90_PUT_ATT(ncid, idlat,'valid_min',MINVAL(rlat) )
  ierr = NF90_PUT_ATT(ncid, idlat,'valid_max',MAXVAL(rlat) )
  ierr = NF90_PUT_ATT(ncid, idlat,'long_name','Latitude')
  ierr = NF90_DEF_VAR(ncid, 'time_counter', NF90_FLOAT, (/ idt /), idtime )
  ierr = NF90_PUT_ATT(ncid, idtime,'units','month since 1900-01-01 00:00:00')
  ierr = NF90_DEF_VAR(ncid, 'socoefr', NF90_FLOAT, (/ idx, idy /)    , idcoef )
  ierr = NF90_PUT_ATT(ncid, idcoef,'valid_min', 0.  )
  ierr = NF90_PUT_ATT(ncid, idcoef,'valid_max', 0.5 )
  ierr = NF90_DEF_VAR(ncid, 'sorunoff', NF90_FLOAT, (/ idx, idy ,idt/), idrnf  )
  ierr = NF90_PUT_ATT(ncid, idrnf,'units','Kg/m2/s')
  ierr = NF90_PUT_ATT(ncid, idrnf,'long_name','runoff flux')
  ierr = NF90_PUT_ATT(ncid, idrnf,'valid_min', 0.  )
  ierr = NF90_PUT_ATT(ncid, idrnf,'valid_min', MAXVAL(sorunoff)  )

  ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL,'Content',' River runoff for main rivers from Dai and Trenberth (2009) ' )
  ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL,'Origin', TRIM(cf_in) )
  ierr = NF90_ENDDEF(ncid)

  ierr = NF90_PUT_VAR(ncid, idlon, rlon )
  ierr = NF90_PUT_VAR(ncid, idlat, rlat )
  ierr = NF90_PUT_VAR(ncid, idtime,  (/ (jm, jm=1,12) /)   )

  ierr = NF90_PUT_VAR(ncid, idcoef, socoefr  )
  ierr = NF90_PUT_VAR(ncid, idrnf , sorunoff )
  ierr = NF90_CLOSE(ncid) 

CONTAINS
  SUBROUTINE GetDaiTrenberth ( ) 
    INTEGER :: js, jt, jy, jmon
    INTEGER :: ijs, ij, iye, iyb, iy, ny
    INTEGER :: ncid, ierr, idt, idstn, id
    INTEGER :: npt
    INTEGER, DIMENSION(:), ALLOCATABLE :: itime
    INTEGER, DIMENSION(12) ::  ncount_month

    REAL(KIND=4) :: rlonmin = -82 , rlonmax = 9.375   ! define window for current model
    REAL(KIND=4) :: rlatmin =  25 , rlatmax = 66      !  
    REAL(KIND=4) :: spval                             ! fill value for flow
    REAL(KIND=4) :: zv                                ! current value

    REAL(KIND=8), DIMENSION(12) :: dsum

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
    ierr = NF90_INQ_VARID (ncid, 'ct_name', id)  ; ierr = NF90_GET_VAR(ncid, id, cname ) ; runoff(:)%ct_name=cname
    ierr = NF90_INQ_VARID (ncid, 'cn_name', id)  ; ierr = NF90_GET_VAR(ncid, id, cname ) ; runoff(:)%cn_name=cname
    ierr = NF90_INQ_VARID (ncid, 'riv_name', id) ; ierr = NF90_GET_VAR(ncid, id, cname ) ; runoff(:)%riv_name=cname
    ierr = NF90_INQ_VARID (ncid, 'ocn_name', id) ; ierr = NF90_GET_VAR(ncid, id, cname ) ; runoff(:)%ocn_name=cname
    DEALLOCATE (cname)


    DO js = 1,npstn
       ALLOCATE ( runoff(js)%flow(npt) )
       ierr = NF90_INQ_VARID( ncid,'station' , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%station , start=(/js/)) !, count=(/1/)   ) 
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 1 ;
       ENDIF
       ierr = NF90_INQ_VARID( ncid,'id'      , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%id      , start=(/js/)) !, count=(/1/)   ) 
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 2 ;
       ENDIF
       ierr = NF90_INQ_VARID( ncid,'lon'     , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%lon     , start=(/js/)) !, count=(/1/)   ) 
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 3 ;
       ENDIF
       ierr = NF90_INQ_VARID( ncid,'lat'     , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%lat     , start=(/js/)) !, count=(/1/)   ) 
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 4 ;
       ENDIF
       ierr = NF90_INQ_VARID( ncid,'lon_mou' , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%lon_mou , start=(/js/)) !, count=(/1/)   ) 
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 5 ;
       ENDIF
       ierr = NF90_INQ_VARID( ncid,'lat_mou' , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%lat_mou , start=(/js/)) !, count=(/1/)   ) 
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 6 ;
       ENDIF
       ierr = NF90_INQ_VARID( ncid,'vol_stn' , id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%vol_stn , start=(/js/)) !, count=(/1/)   ) 
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 7 ;
       ENDIF
       ierr = NF90_INQ_VARID( ncid,'ratio_m2s', id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%ratio_m2s, start=(/js/)) !, count=(/1/) ) 
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 8 ;
       ENDIF
       ierr = NF90_INQ_VARID( ncid,'xnyr',     id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%xnyr    , start=(/js/)) !, count=(/1/)   ) 
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 9 ;
       ENDIF
       ierr = NF90_INQ_VARID( ncid,'yrb',      id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%yrb     , start=(/js/)) !, count=(/1/)   ) 
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 10 ; 
       ENDIF
       ierr = NF90_INQ_VARID( ncid,'yre',      id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%yre     , start=(/js/)) !, count=(/1/)   ) 
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 11 ;
       ENDIF
       ierr = NF90_INQ_VARID( ncid,'elev',     id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%elev    , start=(/js/)) !, count=(/1/)   ) 
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 12 ;
       ENDIF
       ierr = NF90_INQ_VARID( ncid,'FLOW',     id) ; ierr = NF90_GET_VAR(ncid, id, runoff(js)%flow(:),  start = (/ js,1 /), count=(/ 1,npt /) )
       IF ( ierr /= NF90_NOERR )  THEN ; PRINT *, NF90_STRERROR(ierr) ; STOP 18 ;
       ENDIF
       IF ( js == 1 ) ierr = NF90_GET_ATT (ncid, id, '_FillValue', spval )
    ENDDO
    ierr = NF90_CLOSE(ncid)


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
             PRINT *, ' PB with station ', js, TRIM(runoff(js)%riv_name),' ', TRIM(runoff(js)%ct_name), ' ',TRIM(runoff(js)%cn_name), ' ',TRIM(runoff(js)%ocn_name),' month ', jmon, iyb, iye
          ELSE
             runoff(js)%monthly_flow(jmon) = dsum(jmon)/ncount_month(jmon)
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE GetDaiTrenberth
END PROGRAM rnf_compute_area
