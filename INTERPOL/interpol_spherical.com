#!/bin/tcsh -f
##################################################################################################
# Output turbospectrum/babsma format compatible
# A control plot (interpol_check.ps) is displayed at the end.
# Extrapolation is not advised, even if allowed by this program.
# Requires an "cubic" set of 8 MARCS binary format models,
# in other words
# !!!!!   MODELS MUST DIFFER 2 BY 2 BY ONLY ONE PARAMETER !!!!!!
# !!!!!!! ORDER OF THE INPUT MODELS MATTERS !!!!!!!
# here is the order of the files
# model1: Tefflow logglow zlow
# model2: Tefflow logglow zup
# model3: Tefflow loggup zlow
# model4: Tefflow loggup zup
# model5: Teffup logglow zlow
# model6: Teffup logglow zup
# model7: Teffup loggup zlow
# model8: Teffup loggup zup
######################################################################################################

echo 'interpolation script begin'

if (`echo $zref |  awk '{if ($1 < -2.5) print "1";else print "0"}'`) then
set model_path = models/MARCS_st_sph_t02_mod
else
if ($?C) then
# set Csol = `awk '/C /{print $2}' solabu.dat`
# set CFe = `echo $C $Csol $zref | awk '{print $1-$2-$3}'`
# set Cmod =  ( `ls -1d models/c* | xargs -n1 basename | awk '{print ('$CFe'-substr($1,2,40))**2,$1}' | sort -g | awk '{print $2}'` )
set Cmod =  ( `ls -1d models/c* | xargs -n1 basename | awk '{print ('$C'-substr($1,2,40))**2,$1}' | sort -g | awk '{print $2}'` )
set Cmod = $Cmod[1]
else
set Cmod =  'c+0.00'
endif

if ($?alpha) then
set alphamod =   (`ls -1d models/${Cmod}/a* | xargs -n1 basename | awk '{print ('$alpha'-substr($1,2,40))**2,$1}' | sort -g | awk '{print $2}'` )
set alphamod = $alphamod[1]
else
set alphamod = 'a+0.00'
endif

set model_path = models/${Cmod}/${alphamod}/MARCS_st_sph_t02_mod

endif
echo $model_path
#'MARCSbin', 'MARCSweb', 'Atlas' or 'Phoenix' format
set modtype = 'MARCSweb'
#set TEFF = `echo $model | awk '{print substr($1,2,4)}`
#set LOGG  = `echo $model | awk '{print substr($1,8,4)}'`
#set MET =  `echo $model | awk '{print substr($1,26,5)}'`

if (!(-e ${model_path}/grid.list)) then
  ls ${model_path}/s*g*m*z* | xargs -n1 basename | awk '{print $1,substr($1,2,4),substr($1,8,4),substr($1,26,5)}' > ${model_path}/grid.list
  echo grid.list created
endif

source INTERPOL/interpol_spherical_int.com

if ($? == 1 || !(-s ${model}) ) then
    echo 'something went wrong with the interpolation of this model'
    echo ${model}
    echo 'try to extrapolate now'
    source INTERPOL/interpol_ext.com
if ($? == 1 || !(-s ${model}) ) then
exit 1
endif
endif

rm -f models.sm
