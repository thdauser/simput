#! /bin/csh

set simput1 = "athenacrab_test.simput"
set simput2 = "webexample_test.simput"

set all_simput = "$simput1 $simput2"

### 0 ###
# source ${SIMPUT}/bin/simput-install.csh
source ../simput-install.csh

### 1 ### testing simple SIMPUTFILE call
echo "### 1 ### testing simple SIMPUTFILE call ### "
simputfile simput=$simput1 RA=0.0 Dec=0.0 Emin=2. Emax=10. \
    plPhoIndex=2.1 NH=0.4 plFlux=2.052123e-10 clobber=yes


### 2 ### testing simple SIMPUTFILE call
echo "### 1 ### testing WEB EXAMPLE ### "
simputfile RA=40.2 Dec=12.8 XSPECFile="example_spectrum.xcm" \
           LCFile=example_lightcurve.dat MJDREF=50800.0 \
           Emin=0.5 Emax=10.0 srcFlux=2.3e-12 Simput=$simput2

### Cleaning Up
echo " ### Cleaning up ... ###"
if ($1 != "--noclean") then
    rm -vf $all_simput
endif 

echo " ### Test finished! ###"
