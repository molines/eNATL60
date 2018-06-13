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
 1. rnf_xxx tools available in this repository, under TOOLS directory. A Makefile is provided to compile the various fortran programs.
 1. Google Earth highly recommended  in order to acurately locate the rivermouths.

## General procedure:
 * Build a rivermouth file. This file is a 2D file initialized from surface *tmask* in which the zone where to spread the river runoff ('rivermouth') will be  defined, by specifying a unique code number corresponding to each river, identical to the Dai Trenberth coding. 
 * Mask the rivermouth file in order to only keep ocean points.
 * Compute runoff file and check the maximum per cell values ( mm/day). If the run-off exeed 150 mm/day (*arbitrary*), iterate on the rivermouth file in order to increase the spreading of the runoff (increase area).
 

## Roadmap :
### Extending previous rivermouth field (used in NATL60)
 * With NATL60 data file and eNATL60 bathymetry, this is done with the program [data_extend_zone](../TOOLS/data_extend_zone.f90)

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
          variables : Bathymetry !

 * Do it with the rivermouth file used for NATL60. If you do not have this original file, then start from scratch and initialize rivermouth to surface *tmask*

### Prepare kml file with position of Dai-Trenberth river station 
 * This can be done using the bash script [rnf_mk_kml.ksh](../TOOLS/rnf_mk_kml.ksh). The resulting kml file requires some hand editing because of some garbage characters comming from the netcdf Dai and Trenberth netcdf file. 
 * Then in GOOGLE-EARTH you can open this kml file and you will have the Dai-Trenberth station appearing as pin points. This will be very usefull for next step.
 
### Build the rivermouth file using BMGTOOLS.
 * Note that BMGTOOLS requires that the coordinate file has a time axis.
 * Note that BMGTOOLS requires that the variable we are working with is named Bathymetry.
 * For big configuration, it is much easier to work with BMGTOOLS with subdomains. The program [splitfile2](https://github.com/molines/JMMTOOLS/blob/master/TOOLS/splitfile2.f90) can be used to explose the full domain into subdomain and then to merge the subdomain back to full domain. You need to explose both coordinates and data file in the same way :

 >  usage :  splitfile -f IN-file  [-s  x-size y-size] [-n i-size j-size] [-M ] 
       
      PURPOSE :
        This program can be used either for spliting an input file into smaller
        files, whose dimensions are specified, or to merge already splitted 
        files.  The interest of this tool is linked with editing tools  such as
        BMGTOOLS that cannot handle very big files. Each sub-file can then be 
        edited easily and the whole file rebuildt after edition.
            Splitted files have the same root name than the input files, with a
        suffix indicating a rank number, from 1 to the total number of spliited
        files, number 1 corresponding to the south-western most domain, last
        number to the north-eastern most domain.
             When using the merge option, the decomposition is read from the 
        rank 1 file and rebuilt file have IN-file_merged name.
        Split the input files into smaller files with dimensions x-size x y-size
       
      ARGUMENTS :
         -f IN-file : specify the name of the file to split or merge.
                      In case of merge, this is the root name without rank 
                      number. (eg for files like config_0001.nc, just pass 
                      config.)
       
      OPTIONS : 1 and only 1 option must be used.
         -s x-size y-size : specify the size of the subdomain for splitting
         -n i-size j-size : specify the number of subdomains in the I and J 
                            direction.
         -M :  merge splitted files. IN-file is the root name in this case.
         -c coordinates file : use coordinate file to patch nav_lon,nav_lat
              where there are land processors
       
      REQUIRED FILES :
         none
       
      OUTPUT : 
        netcdf files : 
              -s or -n options ROOT_file_nnnn.nc where nnnn is a 4 digit number 
                  indicating the rank of the subdomain. ROOT_file is the name
                  of the input file without the extension.
              -M option : IN_file_merg.nc
           
          variables : same than in IN-file
       
      SEE ALSO :
       BMGTOOLS
 
 * With BMGTOOLS define the zone where to apply the runoff for each river, with the help of google-earth images, to better locate the rivermouth.

 * This procedure may take a lot of time, but make you visit the world !
 
 ### Compute the river runoff
  * This is done with [rnf_compute_runoff](../TOOLS/rnf_compute_runoff.f90) program.
   * in its actual shape, the program assume hard coded name for files :
     - rivermouth.nc
     - coastal-stns-Vol-monthly.updated-Aug2014.nc for data file
     - coordinates.nc 
   * see the code for details and eventual changes !
