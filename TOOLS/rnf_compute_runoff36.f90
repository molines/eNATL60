PROGRAM rnf_compute_runoff36
  !!======================================================================
  !!                     ***  PROGRAM  rnf_compute_runoff36  ***
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
  INTEGER :: narg, iargc, ijarg
  INTEGER :: npi, npj, np, ijs
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
!  CHARACTER(LEN=80) :: cf_in='coastal-stns-Vol-monthly.updated-oct2007.nc'
  CHARACTER(LEN=80) :: cf_in='coastal-stns-Vol-monthly.updated-Aug2014.nc'
  CHARACTER(LEN=80) :: cf_out='runoff.nc'
  CHARACTER(LEN=80) :: cv_e1t = 'e1t'
  CHARACTER(LEN=80) :: cv_e2t = 'e2t'
  CHARACTER(LEN=80) :: cv_gla = 'glamt'
  CHARACTER(LEN=80) :: cv_gph = 'gphit'
  CHARACTER(LEN=80) :: cv_rivmsk = 'Bathymetry'   !  Stupid BMGTOOLS requirement !
  CHARACTER(LEN=80) :: cldum       

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
  LOGICAL            :: lchk = .false.
  !!----------------------------------------------------------------------
  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *, 'USAGE: rnf_compute_runoff36.exe go [-c HGR-file] [-r RIVERMASK-file] '
     PRINT *, '     [-d DATA-runoff-file]'
     PRINT *, ' '
     PRINT *, '   PURPOSE:' 
     PRINT *, '      Compute the runoff file (kg/m2/s) from the river-mask file and'
     PRINT *, '      the Dai/Trenberth data base.'
     PRINT *, ' '
     PRINT *, '   ARGUMENTS:' 
     PRINT *, '      go : if no options are used, any word as argument launch the program.'
     PRINT *, ' '
     PRINT *, '   OPTIONS:' 
     PRINT *, '       [ -c HGR-file] : name of coordinates file for reading the horizontal metrics.'
     PRINT *, '            default is : ',TRIM(cf_coo)
     PRINT *, '       [ -r RIVERMASK-file] : name of the river mask-file'
     PRINT *, '            default is : ',TRIM(cf_riv)
     PRINT *, '       [ -d DATA-runoff-file] : name of Dai/Trenberth data file.' 
     PRINT *, '            default is : ',TRIM(cf_in)
     PRINT *, '    '
     STOP
  ENDIF

  ijarg=1
  DO WHILE (ijarg <= narg )
     CALL getarg( ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE (cldum)
     CASE ( '-c' ) ; CALL getarg(ijarg,cf_coo) ; ijarg=ijarg+1
     CASE ( '-r' ) ; CALL getarg(ijarg,cf_riv) ; ijarg=ijarg+1
     CASE ( '-d' ) ; CALL getarg(ijarg,cf_in ) ; ijarg=ijarg+1
     END SELECT
  ENDDO

  lchk=lchk .OR. chkfile(cf_coo)
  lchk=lchk .OR. chkfile(cf_riv)
  lchk=lchk .OR. chkfile(cf_in )
  IF (lchk) STOP

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
  PRINT *, 'Station  Npts     Area       mm/d(y)  River         mm/d (month)'
  DO js = 1, npstn
     ijs=js
     IF ( js == 1 ) ijs=1000   ! Amazon is tagged as 1000 instead of 1
     np = COUNT ( (itmsk == ijs ) )
     IF ( np /= 0 ) THEN
        itmp(:,:)=0
        WHERE ( itmsk == ijs ) itmp = 1
        darea = SUM( e1t *e2t, (itmp == 1) )
        ! compute runoff ( taking into account some particular case ( St Lawrence, Ottawa, Saguenay for instance)
        SELECT CASE ( ijs )
        CASE ( 1000 ) ! Amazon 
        CASE ( 11)  ! Tocantin + Capim + Guama + Pacajas + Moju
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + &  ! Tocantin
                &    runoff(132)%vol_stn * runoff(132)%ratio_m2s + &   !  Capim
                &    runoff(357)%vol_stn * runoff(357)%ratio_m2s + &   !  Guama
                &    runoff(370)%vol_stn * runoff(370)%ratio_m2s + &   !  Pacajas
                &    runoff(362)%vol_stn * runoff(362)%ratio_m2s       !  Moju
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + &  ! Tocantin
                &            runoff(132)%monthly_flow(:) * runoff(132)%ratio_m2s + &   !  Capim
                &            runoff(357)%monthly_flow(:) * runoff(357)%ratio_m2s + &   !  Guama
                &            runoff(370)%monthly_flow(:) * runoff(370)%ratio_m2s + &   !  Pacajas
                &            runoff(362)%monthly_flow(:) * runoff(362)%ratio_m2s       !  Moju
           runoff( js)%ratio_m2s = 1.  

        CASE ( 16 ) !  ( St Lawrence, Ottawa, Saguenay Ottawa St Maurice Outardes Manicuagan )
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + &
                &    runoff( 59)%vol_stn * runoff( 59)%ratio_m2s + &   !  Ottawa
                &    runoff(130)%vol_stn * runoff(130)%ratio_m2s + &   !  St Maurice
                &    runoff( 66)%vol_stn * runoff( 66)%ratio_m2s    !  Saguenay
!                &    runoff(230)%vol_stn * runoff(230)%ratio_m2s + &   !  Outardes
!                &    runoff(101)%vol_stn * runoff(101)%ratio_m2s       !  Manicuagan
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + &
                &            runoff( 59)%monthly_flow(:) * runoff( 59)%ratio_m2s + &   !  Ottawa
                &            runoff(130)%monthly_flow(:) * runoff(130)%ratio_m2s + &   !  St Maurice
                &            runoff( 66)%monthly_flow(:) * runoff( 66)%ratio_m2s !+ &   !  Saguenay
!                &            runoff(230)%monthly_flow(:) * runoff(230)%ratio_m2s + &   !  Outardes
!                &            runoff(101)%monthly_flow(:) * runoff(101)%ratio_m2s       !  Manicuagan
           runoff( js)%ratio_m2s = 1.  
        CASE ( 29 ) ! Esequibo Cuyuni
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Esequibo
                &    runoff( 93)%vol_stn * runoff( 93)%ratio_m2s                !  Cuyuni
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + &
                &    runoff( 93)%monthly_flow(:) * runoff( 93)%ratio_m2s        !  Cuyuni
           runoff( js)%ratio_m2s = 1.
!        CASE ( 43 ) !  Usumasinta, Grivalva, Rapido de Sama
!           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Usumasinta (MX)
!                &    runoff(133)%vol_stn * runoff(133)%ratio_m2s            + & ! Rapido de Sama
!                &    runoff(201)%vol_stn * runoff(201)%ratio_m2s                ! Grivalva
!           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Usumasinta (MX)
!                &    runoff(133)%monthly_flow(:) * runoff(133)%ratio_m2s    +&  ! Rapido de Sama
!                &    runoff(201)%monthly_flow(:) * runoff(201)%ratio_m2s        ! Grivalva
!           runoff( js)%ratio_m2s = 1.
        CASE ( 45)  ! Rhin + Meuse + Vetch
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + &  ! Rhin
                &    runoff(247)%vol_stn * runoff(247)%ratio_m2s + &   !  Meuse
                &    runoff(606)%vol_stn * runoff(606)%ratio_m2s       !  Vetch
           runoff(js)%monthly_flow = runoff( js)%monthly_flow * runoff( js)%ratio_m2s + &  ! Rhin
                &    runoff(247)%monthly_flow * runoff(247)%ratio_m2s + &   !  Meuse
                &    runoff(606)%monthly_flow * runoff(606)%ratio_m2s       !  Vetch
           runoff( js)%ratio_m2s = 1.  
        CASE ( 47)  ! Caniapiscau + R. aux Melezes
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & !  Caniapiscau
                &    runoff(146)%vol_stn * runoff(146)%ratio_m2s                ! R. aux Melezes
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Caniapiscau
                &    runoff(146)%monthly_flow(:) * runoff(146)%ratio_m2s        ! R. aux Melezes
           runoff( js)%ratio_m2s = 1.
        CASE ( 56)  ! La grande + Sakami
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & !  La grande
                &    runoff(315)%vol_stn * runoff(315)%ratio_m2s                ! Sakami
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! La grande
                &    runoff(315)%monthly_flow(:) * runoff(315)%ratio_m2s        ! Sakami
           runoff( js)%ratio_m2s = 1.
        CASE ( 67)  !  Alabama + Tombigbee
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Alabama
                &    runoff(123)%vol_stn * runoff(123)%ratio_m2s                ! Tombigbee
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & !Alabama
                &    runoff(123)%monthly_flow(:) * runoff(123)%ratio_m2s        ! Tombigbee
           runoff( js)%ratio_m2s = 1.
        CASE ( 74)  ! Susquehanna + Potomac + Rappahannock + James
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + &  ! Susquehanna
                &    runoff(248)%vol_stn * runoff(248)%ratio_m2s + &   !  Potomac
                &    runoff(530)%vol_stn * runoff(530)%ratio_m2s + &   !  Rappahannock
                &    runoff(309)%vol_stn * runoff(309)%ratio_m2s       !  James
           runoff(js)%monthly_flow = runoff( js)%monthly_flow * runoff( js)%ratio_m2s + &  ! Susquehanna
                &    runoff(248)%monthly_flow * runoff(248)%ratio_m2s + &   !  Potomac
                &    runoff(530)%monthly_flow * runoff(530)%ratio_m2s + &   !  Rappahannock
                &    runoff(309)%monthly_flow * runoff(309)%ratio_m2s       !  James
           runoff( js)%ratio_m2s = 1.  
        CASE ( 84)  ! Nottaway + Harricana
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Nottaway
                &    runoff(359)%vol_stn * runoff(359)%ratio_m2s                ! Harricana
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & !Nottaway
                &    runoff(359)%monthly_flow(:) * runoff(359)%ratio_m2s        ! Harricana
           runoff( js)%ratio_m2s = 1.
        CASE ( 88)  ! Garonne + Vezere
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Garonne
                &    runoff(493)%vol_stn * runoff(493)%ratio_m2s                ! Vezere
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Garonne
                &    runoff(493)%monthly_flow(:) * runoff(493)%ratio_m2s        ! Vezere
           runoff( js)%ratio_m2s = 1.
        CASE (114)  ! Douro + Tamega
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Douro
                &    runoff(463)%vol_stn * runoff(463)%ratio_m2s                ! Tamega
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Douro
                &    runoff(463)%monthly_flow(:) * runoff(463)%ratio_m2s        ! Tamega
           runoff( js)%ratio_m2s = 1.

        CASE (116)  ! Kouilou +Niari
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Kouilou
                &    runoff(235)%vol_stn * runoff(235)%ratio_m2s                ! Niari
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & !Kouilou
                &    runoff(235)%monthly_flow(:) * runoff(235)%ratio_m2s        ! Niari
           runoff( js)%ratio_m2s = 1.
        CASE (132 ) ! Guama + Capim
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Guama
                &    runoff(357)%vol_stn * runoff(357)%ratio_m2s                ! Capim
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & !Guama
                &    runoff(357)%monthly_flow(:) * runoff(357)%ratio_m2s        ! Capim
           runoff( js)%ratio_m2s = 1.
        CASE (133)  ! Rompido de Sama + Grijalva
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Guama
                &    runoff(201)%vol_stn * runoff(201)%ratio_m2s                ! Capim
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & !Guama
                &    runoff(201)%monthly_flow(:) * runoff(201)%ratio_m2s        ! Capim
           runoff( js)%ratio_m2s = 1.
        CASE (137)  !  Glomma + Dramselv
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Glomma
                &    runoff(256)%vol_stn * runoff(256)%ratio_m2s                ! Dramselv
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Glomma
                &    runoff(256)%monthly_flow(:) * runoff(256)%ratio_m2s        ! Dramselv
           runoff( js)%ratio_m2s = 1.
        CASE (141)  ! Daugava + Lielupe + Aiviekste
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + &  ! Daugava
                &    runoff(500)%vol_stn * runoff(500)%ratio_m2s + &   !  Lielupe
                &    runoff(540)%vol_stn * runoff(540)%ratio_m2s       !  Aiviekste
           runoff(js)%monthly_flow = runoff( js)%monthly_flow * runoff( js)%ratio_m2s + &  ! Daugava
                &    runoff(500)%monthly_flow * runoff(500)%ratio_m2s + &   !  Lielupe
                &    runoff(540)%monthly_flow * runoff(540)%ratio_m2s       !  Aiviekste
           runoff( js)%ratio_m2s = 1.  
        CASE (170)  ! Itapecuru + Munim
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Itapecuru
                &    runoff(374)%vol_stn * runoff(374)%ratio_m2s                ! Munim
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Itapecuru
                &    runoff(374)%monthly_flow(:) * runoff(374)%ratio_m2s        ! Munim
           runoff( js)%ratio_m2s = 1.
        CASE (179)  ! Brazos + Trinity
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Brazos
                &    runoff(295)%vol_stn * runoff(295)%ratio_m2s                ! Trinity
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Brazos
                &    runoff(295)%monthly_flow(:) * runoff(295)%ratio_m2s        ! Trinity
           runoff( js)%ratio_m2s = 1.
        CASE (181)  ! Jaguaribe + Bababuiu
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Jaguaribe
                &    runoff(520)%vol_stn * runoff(520)%ratio_m2s                ! Bababuiu
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Jaguaribe
                &    runoff(520)%monthly_flow(:) * runoff(520)%ratio_m2s        ! Bababuiu
           runoff( js)%ratio_m2s = 1.
        CASE (227)  ! Hvita + Bruara
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Hvita
                &    runoff(583)%vol_stn * runoff(583)%ratio_m2s                ! Bruara
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Hvita
                &    runoff(583)%monthly_flow(:) * runoff(583)%ratio_m2s        ! Bruara
           runoff( js)%ratio_m2s = 1.
        CASE (288)  ! Sabines + Noches
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Sabines
                &    runoff(317)%vol_stn * runoff(317)%ratio_m2s                ! Noches
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Sabines
                &    runoff(317)%monthly_flow(:) * runoff(317)%ratio_m2s        ! Noches
           runoff( js)%ratio_m2s = 1.
        CASE (306)  !Pindare + Mearim + Grajau + Turiacu
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Pindare
                &    runoff(395)%vol_stn * runoff(395)%ratio_m2s            + & ! Mearim
                &    runoff(397)%vol_stn * runoff(397)%ratio_m2s            + & ! Grajau
                &    runoff(445)%vol_stn * runoff(445)%ratio_m2s                ! Turiacu
           runoff(js)%monthly_flow = runoff( js)%monthly_flow * runoff( js)%ratio_m2s + & ! Pindare
                &    runoff(395)%monthly_flow * runoff(395)%ratio_m2s        + & ! Mearim
                &    runoff(397)%monthly_flow * runoff(397)%ratio_m2s        + & ! Grajau
                &    runoff(445)%monthly_flow * runoff(445)%ratio_m2s            ! Turiacu
           runoff( js)%ratio_m2s = 1.
        CASE (417)  ! Nidelv + Gaula
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Nidelv
                &    runoff(462)%vol_stn * runoff(462)%ratio_m2s                ! Gaula
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Nidelv
                &    runoff(462)%monthly_flow(:) * runoff(462)%ratio_m2s        ! Gaula
           runoff( js)%ratio_m2s = 1.
        CASE (423 ) ! Ouergha + Sehou
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Ouergah
                &    runoff(498)%vol_stn * runoff(498)%ratio_m2s                ! Sehou
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Ouergha
                &    runoff(498)%monthly_flow(:) * runoff(498)%ratio_m2s        ! Sehou
           runoff( js)%ratio_m2s = 1.
        CASE (458) ! Wye + severn
           runoff(js)%vol_stn = runoff( js)%vol_stn * runoff( js)%ratio_m2s + & ! Wye
                &    runoff(486)%vol_stn * runoff(486)%ratio_m2s                ! Severn
           runoff(js)%monthly_flow(:) = runoff( js)%monthly_flow(:) * runoff( js)%ratio_m2s + & ! Wye
                &    runoff(486)%monthly_flow(:) * runoff(486)%ratio_m2s        ! Severn
           runoff( js)%ratio_m2s = 1.
        END SELECT
!        PRINT *, js, np , darea/1.e6, 'km^2', runoff(js)%vol_stn * runoff(js)%ratio_m2s/darea * dconvcoef*dconv2mm, TRIM(runoff(js)%riv_name), SUM(runoff(js)%monthly_flow(:))/12.*runoff(js)%ratio_m2s/darea*1000.*86400. ! mm/day
        PRINT 777 ,  js, np , darea/1.e6, 'km^2', runoff(js)%vol_stn * runoff(js)%ratio_m2s/darea * dconvcoef*dconv2mm, &
                     TRIM(runoff(js)%riv_name), SUM(runoff(js)%monthly_flow(:))/12.*runoff(js)%ratio_m2s/darea*1000.*86400. ! mm/day
777 FORMAT(i5,i8, f10.2,x,a,f8.2,x,a, f8.2)

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
  LOGICAL FUNCTION chkfile (cd_file)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION chkfile  ***
    !!
    !! ** Purpose :  Check if cd_file exists.
    !!               Return false if it exists, true if it does not
    !!               Do nothing is filename is 'none'
    !!
    !! ** Method  : Doing it this way allow statements such as
    !!              IF ( chkfile( cf_toto) ) STOP 99  ! missing file
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),  INTENT(in) :: cd_file

    INTEGER(KIND=4)               :: ierr
    LOGICAL                       :: ll_exist
    !!----------------------------------------------------------------------
    IF ( TRIM(cd_file) /= 'none')  THEN
       INQUIRE (file = TRIM(cd_file), EXIST=ll_exist)

       IF (ll_exist) THEN
          chkfile = .false.
       ELSE
          PRINT *, ' File ',TRIM(cd_file),' is missing '
          chkfile = .true.
       ENDIF
    ELSE
       chkfile = .false.  ! 'none' file is not checked
    ENDIF

  END FUNCTION chkfile

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
END PROGRAM rnf_compute_runoff36
