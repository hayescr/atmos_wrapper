
set Teff_list2 = (`awk '{print $0}' ${model_path}/grid.list | awk '{print $2}' | sort -u`)
set Teff_list = (`echo $Teff_list2 | awk 'BEGIN{RS=" "}{printf "%f %s \n",(('${Tref}' - $1))**2,$1}'  | sort -u | sort -n | awk '{if ($1 == 0) print $2,$2; else if (NR == 1 || NR == 2) print $2}'`)
set Teff_min = `echo $Teff_list | awk '{if ($1 > $2) print $2;else print $1 }' `
set Teff_max =  `echo $Teff_list | awk '{if ($1 < $2) print $2;else print $1 }'`
echo "Teff min" $Teff_min "max" $Teff_max
##
set logg_list1 = (`awk '{if ($2 == '$Teff_min') print $0}' ${model_path}/grid.list | awk '{print $3}' | sort -u`)
set logg_list2 = ()
foreach loggmetmin ($logg_list1)
set logg_list2 = (`awk '{if ($2 == '$Teff_max' && $3 == '$loggmetmin') print $0}' ${model_path}/grid.list | awk '{print $3}' | sort -u` $logg_list2)
end
set logg_list = (`echo $logg_list2 | awk 'BEGIN{RS=" "}{printf "%f %s \n",(('${loggref}' - $1))**2,$1}'  | sort -u | sort -n | awk '{if ($1 == 0) print $2,$2; else if (NR == 1 || NR == 2) print $2}'`)
set logg_min = `echo $logg_list | awk '{if ($1 > $2) print $2;else print $1 }' `
set logg_max =  `echo $logg_list | awk '{if ($1 < $2) print $2;else print $1 }'`
echo "logg min" $logg_min "max" $logg_max
##
set met_list1 =  (`awk '{if ($2 == '$Teff_min' && $3 == '$logg_min') print $4}' ${model_path}/grid.list | sort -u`)
set met_list2 = ()
foreach metteffmin ($met_list1)
set met_list2 =  (`awk '{if ($2 == '$Teff_max' && $3 == '$logg_min' && $4 == '$metteffmin') print $4}' ${model_path}/grid.list | sort -u` $met_list2)
end
set met_list3 = ()
foreach metteffmin ($met_list2)
set met_list3 =  (`awk '{if ($2 == '$Teff_min' && $3 == '$logg_max' && $4 == '$metteffmin') print $4}' ${model_path}/grid.list | sort -u` $met_list3)
end
set met_list4 = ()
foreach metteffmin ($met_list3)
set met_list4 =  (`awk '{if ($2 == '$Teff_max' && $3 == '$logg_max' && $4 == '$metteffmin') print $4}' ${model_path}/grid.list | sort -u` $met_list4)
end
set met_min = (` echo $met_list4 | awk 'BEGIN{RS=" "}{if (('${zref}' - $1) >=0) print ('${zref}' - $1),$1}' | sort -u | sort -n | awk '{if (NR == 1) print $2}' `)
if ($#met_min == 0) then
set met_min = (`  echo $met_list4 | awk 'BEGIN{RS=" "}{if (('${zref}' - $1) <0) print ('${zref}' - $1)**2,$1}' | sort -u | sort -n | awk '{if (NR == 1) print $2}' `)
set met_max = (`  echo $met_list4 | awk 'BEGIN{RS=" "}{if (('${zref}' - $1) <0) print ('${zref}' - $1)**2,$1}' | sort -u | sort -n | awk '{if (NR == 2) print $2}' `)
else
set met_max = (`  echo $met_list4 | awk 'BEGIN{RS=" "}{if (('${zref}' - $1) <0) print ($1 - '${zref}'),$1}' | sort -u | sort -n | awk '{if (NR == 1) print $2}'`)
if ($#met_max == 0) then
set met_min = (`  echo $met_list4 | awk 'BEGIN{RS=" "}{if (('${zref}' - $1) >= 0) print ($1 - '${zref}')**2,$1}' | sort -u | sort -n | awk '{if (NR == 2) print $2}'`)
set met_max = (` echo $met_list4 | awk 'BEGIN{RS=" "}{if (('${zref}' - $1) >= 0) print ($1 - '${zref}')**2,$1}' | sort -u | sort -n | awk '{if (NR == 1) print $2}'`)
endif
endif
if (`echo $met_min  $zref | awk '{if ($1 == $2) print "1";else print "0"}'`) then
set met_max = $met_min
endif
if (`echo $met_max  $zref | awk '{if ($1 == $2) print "1";else print "0"}'`) then
set met_min = $met_max
endif
echo "z min" $met_min "max" $met_max


######
set model1 = (`awk '{if ($2 == '$Teff_min' && $3 == '$logg_min' && $4 == '$met_min') print $1}'  ${model_path}/grid.list`)
set model1 = $model1[$#model1]
 echo $model1
if (${model1} == "") then
echo model1 does not exist
set failed = 1
goto failed
endif

set model2 = (`awk '{if ($2 == '$Teff_min' && $3 == '$logg_min' && $4 == '$met_max' ) print $1}'  ${model_path}/grid.list`)
set model2 = $model2[$#model2]
 echo $model2
if (${model2} == "") then
echo model2 does not exist
set failed = 1
goto failed
endif

set model3 = (`awk '{if ($2 == '$Teff_min' && $3 == '$logg_max' && $4 == '$met_min') print $1}'  ${model_path}/grid.list`)
set model3 = $model3[$#model3]
 echo $model3
if (${model3} == "") then
echo model3 does not exist
set failed = 1
goto failed
endif

set model4 = (`awk '{if ($2 == '$Teff_min' && $3 == '$logg_max' && $4 == '$met_max' ) print $1}'  ${model_path}/grid.list`)
set model4 = $model4[$#model4]
 echo $model4
if (${model4} == "") then
echo model4 does not exist
set failed = 1
goto failed
endif


set model5 = (`awk '{if ($2 == '$Teff_max' && $3 == '$logg_min' && $4 == '$met_min'  ) print $1}'  ${model_path}/grid.list`)
set model5 = $model5[$#model5]
 echo $model5
if (${model5} == "") then
echo model5 does not exist
set failed = 1
goto failed
endif

set model6 = (`awk '{if ($2 == '$Teff_max' && $3 == '$logg_min' && $4 == '$met_max' ) print $1}'  ${model_path}/grid.list`)
set model6 = $model6[$#model6]
 echo $model5
if (${model6} == "") then
echo model6 does not exist
set failed = 1
goto failed
endif

set model7 = (`awk '{if ($2 == '$Teff_max' && $3 == '$logg_max' && $4 == '$met_min'  ) print $1}'  ${model_path}/grid.list`)
set model7 = $model7[$#model7]
 echo $model7
if (${model7} == "") then
echo model7 does not exist
set failed = 1
goto failed
endif

set model8 = (`awk '{if ($2 == '$Teff_max' && $3 == '$logg_max' && $4 == '$met_max' ) print $1}'  ${model_path}/grid.list`)
set model8 = $model8[$#model8]
 echo $model8
if (${model8} == "") then
echo model8 does not exist
set failed = 1
goto failed
endif


#enter here the values requested for the interpolated model
set model =  ${rootmodel}.int
echo $model

#format of output: 'TurboSpectrum',  'Atlas', 'MOOG',default
if ($?outtype) then
  echo $outtype
else
  set outtype = 'TurboSpectrum'
endif
#### the test option is set to .true. if you want to plot comparison model (model_test)
set test = '.false.'
set model_test = ' /home/tpm40/data/MARCS/models/CS29512-073/5500g3.2z-2.0a0.4CN'
#'MARCSbin', 'MARCSweb', 'Atlas' or 'Phoenix' format
set model_test_type = 'MARCSbin'
# interpolation program (for further details see interpol_models.f)
INTERPOL/interpol_modeles <<EOF
'${model_path}/${model1}'
'${model_path}/${model2}'
'${model_path}/${model3}'
'${model_path}/${model4}'
'${model_path}/${model5}'
'${model_path}/${model6}'
'${model_path}/${model7}'
'${model_path}/${model8}'
'${model}'
${modtype}
${Tref}
${loggref}
${zref}
${outtype}
${test}
$model_test
$model_test_type
EOF
if (-e ${rootmodel}.ext) then
set model = ${rootmodel}.ext
endif

failed:
