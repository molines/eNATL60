# Preparing runoff file for eNATL60
  In this simulation, we still emulate the river runoff as a precipitation spread over the rivermouth. This procedure requires to have a 2D file with 2 variables: **sorunoff**  giving the amount of 'precip' (kg/m2/s) to add at the rivermouth. This sorunoff variable is in general a monthly climatology. **socoefr** representing the rivermouthes and used as a mask field for specific actions such as for example, shutting off the SSS restoring on the rivermouth. For historical reason this socoefr variable has values between 0 and 0.5 (In fact either 0 or 0.5).

In this document we present the procedure used for producing the runoff file.

## Required information :
 1. Coordinate file for the domain
 1. Surface mask file.
 1. Bathymetry file of the domain (finaly used as a mask)
 1. [Dai and Trenberth dataset](http://www.cgd.ucar.edu/cas/catalog/surface/dai-runoff/)
   * we use the last update of the monthly dataset (coastal-stns-Vol-monthly.updated-Aug2014.nc)

## Required software :
 1. [BIMGTOOLS](http://archimer.ifremer.fr/doc/00195/30646/) for editing the rivermouth file
 1. rnf_xxx tools available in this repository.
 1. Google Earth highly recommended  in order to acurately locate the rivermouths.

## General procedure:
 * Build a rivermouth file. This file is a 2D file initialized from surface *tmask* in which the zone where to spread the river runoff ('rivermouth') will be  defined, by specifying a unique code number corresponding to each river, identical to the Dai Trenberth coding. 
 * Mask the rivermouth file in order to only keep ocean points.
 * Compute runoff file and check the maximum per cell values ( mm/day). If the run-off exeed 150 mm/day (*arbitrary*), iterate on the rivermouth file in order to increase the spreading of the runoff (increase area).
 

## Roadmap :
### Extending previous rivermouth field (used in NATL60)
 * With NATL60 data file and eNATL60 bathymetry, this is done with the program [data_extend_zone](TOOLS/data_extend_zone.f90)

> usage :  data_extend_zone -b BATHY-new -d DATA-old -v VAR-old ...
           [-o FILE-out] [-w imin imax jmin jmax ]
     
      PURPOSE :
         Create a file of the size of the new bathymetry, with old values copied
          from old data file
      
      ARGUMENTS :
        -b BATHY-new : give the name of the new bathymetry 
        -d DATA-old  : give the name of the old data file
        -v VAR-old  : give the name of the variable to copy in the old data file
     
      OPTIONS :
        -w imin imax jmin jmax : specify the position of old grid with respect 
           to the new one. Default is for NATL60 into eNATL60
        -o FILE-out: specify name of output file, instead of <BATHY-new>.wrk   
    
      OUTPUT : 
        netcdf file : none
          variables : Bathymetry !!

 * Do it with the rivermouth file used for NATL60. If you do not have this original file, then start from scratch and initialize rivermouth to surface *tmask*

### Prepare kml file with position of Dai-Trenberth river station 


