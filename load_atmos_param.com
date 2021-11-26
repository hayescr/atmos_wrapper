#!/bin/tcsh -f

set debug = ''
set debug = '1'
if ($?debug) then
  if ($debug) then
    set logfile = /dev/stdout
  else
    set logfile = /dev/null
  endif
else
set logfile = /dev/stdout
endif

set Tref = $1
set loggref = $2
set zref = $3
set TURBVEL = $4
set C = $5
set alpha = $6
set outtype = $7

if ("$8" == "") then
echo "No filename provided"
else
set model = $8
endif

echo "stellar parameters selected (Teff,logg,[Fe/H],microturb):" $Tref $loggref $zref $TURBVEL


if (`echo $loggref | awk '{if ($1 < 3) print "1";else print "0"}'`) then
  # set model = ${star}.int
  # set model = ${Tref}g${loggref}m1.0z${zref}.int
  set SPH = T
  source INTERPOL/interpol_spherical.com  >& $logfile
else
  # set model = ${Tref}g${loggref}z${zref}_${star}.int
  # set model = ${Tref}g${loggref}z${zref}.int
  set SPH = F
  source INTERPOL/interpol_planparallel.com  >& $logfile
endif

if ($? == 1  || !(-s ${model})) then
  echo 'something went wrong with the interpolation of this model'
  echo ${model}
  exit 1
endif
