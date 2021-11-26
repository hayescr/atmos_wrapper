echo ${model_path}
set met_min = (` awk '{if (('${zref}' - $4) >=0) print ('${zref}' - $4),$4}' ${model_path}/grid.list | sort -u | sort -n | awk '{if (NR == 1) print $2}' `)
if ($#met_min == 0) then
set met_min = (` awk '{if (('${zref}' - $4) <0) print ('${zref}' - $4)**2,$4}' ${model_path}/grid.list | sort -u | sort -n | awk '{if (NR == 1) print $2}' `)
set met_max = (` awk '{if (('${zref}' - $4) <0) print ('${zref}' - $4)**2,$4}' ${model_path}/grid.list | sort -u | sort -n | awk '{if (NR == 2) print $2}' `)
else
set met_max = (` awk '{if (('${zref}' - $4) <0) print ($4 - '${zref}'),$4}' ${model_path}/grid.list | sort -u | sort -n | awk '{if (NR == 1) print $2}'`)
if ($#met_max == 0) then
set met_min = (` awk '{if (('${zref}' - $4) >= 0) print ($4 - '${zref}')**2,$4}' ${model_path}/grid.list | sort -u | sort -n | awk '{if (NR == 2) print $2}'`)
set met_max = (` awk '{if (('${zref}' - $4) >= 0) print ($4 - '${zref}')**2,$4}' ${model_path}/grid.list | sort -u | sort -n | awk '{if (NR == 1) print $2}'`)
endif
endif
if (`echo $met_min  $zref | awk '{if ($1 == $2) print "1";else print "0"}'`) then
set met_max = $met_min
endif
if (`echo $met_max  $zref | awk '{if ($1 == $2) print "1";else print "0"}'`) then
set met_min = $met_max
endif
echo "z min" $met_min "max" $met_max "zref=" $zref

set Teff_list1 = (`awk '{if ($4 == '$met_min') print $0}' ${model_path}/grid.list | awk '{print $2}' | sort -u`)
set Teff_list2 = ()
foreach Tmetmin ($Teff_list1)
set Teff_list2 = (`awk '{if ($4 == '$met_max' && $2 == '$Tmetmin') print $0}' ${model_path}/grid.list | awk '{print $2}' | sort -u` $Teff_list2)
end
set Teff_list = (`echo $Teff_list2 | awk 'BEGIN{RS=" "}{printf "%f %s \n",(('${Tref}' - $1))**2,$1}'  | sort -u | sort -n | awk '{if ($1 == 0) print $2,$2; else if (NR == 1 || NR == 2) print $2}'`)
set Teff_min = `echo $Teff_list | awk '{if ($1 > $2) print $2;else print $1 }' `
set Teff_max =  `echo $Teff_list | awk '{if ($1 < $2) print $2;else print $1 }'`
echo "Teff min" $Teff_min "max" $Teff_max

set logg_list1 = (`awk '{if ($2 == '$Teff_min' && $4 == '$met_min') print $3}' ${model_path}/grid.list`)
set logg_list2 = ()
foreach loggmetmin ($logg_list1)
set logg_list2 = (`awk '{if ($2 == '$Teff_max' && $4 == '$met_min' && $3 == '$loggmetmin') print $3}' ${model_path}/grid.list | sort -u` $logg_list2)
end
set logg_list3 = ()
foreach loggmetmin ($logg_list2 )
set logg_list3 = (`awk '{if ($2 == '$Teff_min' && $4 == '$met_max' && $3 == '$loggmetmin') print $3}' ${model_path}/grid.list | sort -u` $logg_list3)
end
set logg_list4 = ()
foreach loggmetmin ($logg_list3 )
set logg_list4 = (`awk '{if ($2 == '$Teff_max' && $4 == '$met_max' && $3 == '$loggmetmin') print $3}' ${model_path}/grid.list | sort -u` $logg_list4)
end
set logg_list = (`echo  $logg_list4 | awk 'BEGIN{RS=" "}{print ('${loggref}' - $1)**2,$1}' | sort -u | sort -n | awk '{if ($1 == 0) print $2,$2; else if (NR == 1 || NR == 2) print $2}'`)
set logg_min = `echo $logg_list | awk '{if ($1 > $2) print $2;else print $1 }' `
set logg_max =  `echo $logg_list | awk '{if ($1 < $2) print $2;else print $1 }'`
echo "logg min" $logg_min "max" $logg_max

set model1 = (`awk '{if ($2 == '$Teff_min' && $3 == '$logg_min' && $4 == '$met_min') print $1}' ${model_path}/grid.list `)
set model1 = $model1[$#model1]
 echo $model1
if (${model1} == "") then
echo model1 does not exist
set failed = 1
goto failed
endif

set model2 = (`awk '{if ($2 == '$Teff_min' && $3 == '$logg_min' && $4 == '$met_max' ) print $1}' ${model_path}/grid.list `)
set model2 = $model2[$#model2]
 echo $model2
if (${model2} == "") then
echo model2 does not exist
set failed = 1
goto failed
endif

set model3 = (`awk '{if ($2 == '$Teff_min' && $3 == '$logg_max' && $4 == '$met_min') print $1}' ${model_path}/grid.list `)
set model3 = $model3[$#model3]
 echo $model3
if (${model3} == "") then
echo model3 does not exist
set failed = 1
goto failed
endif

set model4 = (`awk '{if ($2 == '$Teff_min' && $3 == '$logg_max' && $4 == '$met_max' ) print $1}' ${model_path}/grid.list `)
set model4 = $model4[$#model4]
 echo $model4
if (${model4} == "") then
echo model4 does not exist
set failed = 1
goto failed
endif


set model5 = (`awk '{if ($2 == '$Teff_max' && $3 == '$logg_min' && $4 == '$met_min'  ) print $1}' ${model_path}/grid.list `)
set model5 = $model5[$#model5]
 echo $model5
if (${model5} == "") then
echo model5 does not exist
set failed = 1
goto failed
endif

set model6 = (`awk '{if ($2 == '$Teff_max' && $3 == '$logg_min' && $4 == '$met_max' ) print $1}' ${model_path}/grid.list `)
set model6 = $model6[$#model6]
 echo $model5
if (${model6} == "") then
echo model6 does not exist
set failed = 1
goto failed
endif

set model7 = (`awk '{if ($2 == '$Teff_max' && $3 == '$logg_max' && $4 == '$met_min'  ) print $1}' ${model_path}/grid.list `)
set model7 = $model7[$#model7]
 echo $model7
if (${model7} == "") then
echo model7 does not exist
set failed = 1
goto failed
endif

set model8 = (`awk '{if ($2 == '$Teff_max' && $3 == '$logg_max' && $4 == '$met_max' ) print $1}' ${model_path}/grid.list `)
set model8 = $model8[$#model8]
 echo $model8
if (${model8} == "") then
echo model8 does not exist
set failed = 1
goto failed
endif



#enter here the values requested for the interpolated model
set rootmodel = `echo ${model} | awk '{print substr($1,1,index($1,".int")-1)}'`
set failed = 0
set model_done = 1
echo $model

#format of output: 'TurboSpectrum',  'Atlas', 'MOOG',default
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
if (!($failed)) then
rm -f ${model}_record models_record.sm
cp -f $model ${model}_record
cp -vf models.sm models_record.sm
set initmodel_record = ${model}
set failed_previous = 0
else
set failed_previous = 1
endif

if ($failed || `echo $model | awk '{if ($1 ~ "ext") print "1";else print "0"}'` ) then
echo "extrapolated model; try alternate3"
rm -f ${model}
set failed = 0

source INTERPOL/interpol_alternate2.com
echo alternate2 done echo $model
else
source INTERPOL/interpol_alternate3.com
echo alternate3 done echo $model
endif
if (!($failed) && $failed_previous) then
rm -f ${model}_record  models_record.sm
cp -f $model ${model}_record
cp -vf models.sm models_record.sm
set initmodel_record = ${model}
set failed_previous = 0
else
if ($failed && $failed_previous) then
set failed_previous = 1
endif

if ($failed || `echo $model | awk '{if ($1 ~ "ext") print "1";else print "0"}'` ) then
echo "still extrapolated model; try alternate1"
rm -f ${model}
set failed = 0
source INTERPOL/interpol_alternate1.com
if (!($failed) && $failed_previous) then
rm -f ${model}_record models_record.sm
cp -f $model ${model}_record
cp -vf models.sm models_record.sm
set initmodel_record = ${modele}
set failed_previous = 0
else
if ($failed && $failed_previous) then
set failed_previous = 1
endif
endif
if ($failed || `echo $model | awk '{if ($1 ~ "ext") print "1";else print "0"}'` ) then
echo "still extrapolated model; try last alternate"
rm -f ${model}
set failed = 0
source INTERPOL/interpol_alternate3.com
else
source INTERPOL/interpol_alternate2.com
endif
if (!($failed) && $failed_previous) then
rm -f ${model}_record models_record.sm
cp -f $model ${model}_record
cp -vf models.sm models_record.sm
set initmodel_record = ${model}
set failed_previous = 0
else
if ($failed && $failed_previous) then
set failed_previous = 1
endif


endif
endif
endif

if (($failed || `echo $model | awk '{if ($1 ~ "ext") print "1";else print "0"}'` ) && ! ($failed_previous)) then
rm -f ${model}
echo "no alternate worked; back to initial extrapolated"
set model =  $initmodel_record
rm -f models.sm
cp -f ${model}_record  $model
cp -vf models_record.sm models.sm
else
exit 1
endif
skipmodel:


rm -f  models_record.sm  ${model}_record
