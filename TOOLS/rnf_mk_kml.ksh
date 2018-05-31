#!/bin/ksh
# This script build a kml file (suitable for google earth) with push-pins located at  Dai_Trenberth stations.


# rnf_reader read the data base and output essencial fields (tab separated).
./rnf_reader.exe > ztmp_rivers.txt

cat rnf_header.tmpl > Dai_TrenBerth.kml
# there are a lot of garbage characters that are filtered with the gsub.
cat ztmp_rivers.txt | awk  -F '\t' '{ if (NR > 1 ) {rid=$1 ; lon=$3 ; lat = $4 ; vol=$5  ; riv = $9   ; stn = $10 ; \
      gsub(/\0/,"",riv);  \
      gsub(/[\000\370\300\360\231\024\314@\212\200\033\240\340\363[:punct:]\304\315\320\262\256\321\337\260\237\350\246\354\003\346\341\263]/,"",stn) ; \
      system ( "sed -e \"s/LON/"lon"/g\"  -e \"s/LAT/"lat"/g\"   -e \"s/COMMENT/"rid"     "vol "km3 per yr   " stn"/g\"  -e \"s/NAME/"riv"/g\" rnf_placemark.tmpl  >> Dai_TrenBerth.kml" ) } }' 

cat rnf_coda.tmpl >> Dai_TrenBerth.kml
