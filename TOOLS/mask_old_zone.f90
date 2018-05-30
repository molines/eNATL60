PROGRAM mask_old_zone
  !!======================================================================
  !!                     ***  PROGRAM  mask_old_zone  ***
  !!=====================================================================
  !!  ** Purpose : modify the coast-line of a newly produced bathymetry
  !!               by create_bathy.exe (Nesting Tools), by applying a mask 
  !!               of an already existing bathymetry in a subdomain.
  !!           Example: set the NATL60 coast-line in eNATL60.
  !!
  !!  ** Method  : Position of subdomain is hard coded via the coordinates of
  !!               the corners. All sea points of the old bathymetry are kept.
  !!               In case of a dry new point correponding to a wet old point,
  !!               bathymetry is set to rdepmin ( hardcoded to 3m here).
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
  INTEGER, DIMENSION(4) :: icorner, jcorner
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: mask        ! on eDomain, inner part on Domain

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: bathy  ! on eDomain
  REAL(KIND=4)                              :: rdepmin=3 !m

  CHARACTER(LEN=80) :: cf_bathy='eNATL60_BATHY_GEBCO_2014_2D_raw_v2.nc'
  CHARACTER(LEN=80) :: cf_mask='NATL60_v4.1_cdf_byte_mask.nc'

! initialise corner value : ie index of Domain corners in eDomain reference

!      4 --------  3
!      |           |
!      |           |
!      1 --------  2
!  icorner=(/1053,6474,6474,1053/)   ! v0-v1
!  jcorner=(/ 641, 641,4094,4094/)   ! v0-v1
  icorner=(/1053,6474,6474,1053/)    ! v2
  jcorner=(/1276,1276,4729,4729/)    ! v2

  nx=icorner(2) - icorner (1) + 1
  ny=jcorner(4) - jcorner (1) + 1
   

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
  ierr=NF90_GET_VAR(ncidm, id, mask(icorner(1):icorner(2), jcorner(1):jcorner(4)), start=(/1,1,1,1/),count=(/nx,ny,1,1/) )
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

END PROGRAM mask_old_zone

