

set met_min = (`cat ${model_path}/grid.list  | awk '{print substr($1,index($1,"z")+1,40)}' | awk '{print substr($1,1,index($1,"a")-2)}' | sort -u | sort -n | awk '{if (('${zref}' - $1) >=0) print ('${zref}' - $1),$1}' | sort -n | awk '{if (NR == 1) print $2}' `)
set met_max = (`cat ${model_path}/grid.list | awk '{print substr($1,index($1,"z")+1,40)}' | awk '{print substr($1,1,index($1,"a")-2)}' | sort -u | sort -n | awk '{if (('${zref}' - $1) <=0) print ($1 - '${zref}'),$1}' | sort -n | awk '{if (NR == 1) print $2}' `)

if ($#met_min == 0 || $#met_max == 0) then
  echo extrapolation
  exit 1
endif

echo met = $zref $met_min $met_max

set logg_min = (`awk '/p*.*g.*m.*.0.*z\'$met_min'.*/{print $0}' ${model_path}/grid.list | awk '{print substr($1,8,index($1,"m")-9)}' | sort -u | sort -n | awk '{if (('${loggref}' - $1) >= 0) print ('${loggref}' - $1),$1}' | sort -n | awk '{if (NR <2) print $2}'`)
set logg_max = (`awk '/p*.*g.*m.*.0.*z\'$met_min'.*/{print $0}' ${model_path}/grid.list  | awk '{print substr($1,8,index($1,"m")-9)}' | sort -u | sort -n | awk '{if (('${loggref}' - $1) <= 0) print ($1- '${loggref}'),$1}' | sort -n | awk '{if (NR <2) print $2}'`)

if ($#logg_min == 0 || $#logg_max == 0) then
  echo extrapolation
  exit 1
endif

echo logg = $loggref $logg_min $logg_max


set Teff_min = (`awk '/p.*g\'$logg_min'.*m.*.0.*z\'$met_min'.*/{print $0}' ${model_path}/grid.list | awk '{print substr($1,2,4)}' | sort -u | sort -n| awk '{if (('${Tref}' - $1) >= 0) print ('${Tref}' - $1),$1}' | sort -n | awk '{if (NR <2) print $2}' `)

set Teff_max = (`awk '/p.*g\'$logg_min'.*m.*.0.*z\'$met_min'.*/{print $0}' ${model_path}/grid.list | awk '{print substr($1,2,4)}' | sort -u | sort -n| awk '{if (('${Tref}' - $1) <= 0) print ($1 - '${Tref}'),$1}' | sort -n | awk '{if (NR <2) print $2}' `)

if ($#Teff_min == 0 || $#Teff_max == 0) then
  echo extrapolation
  exit 1
endif

echo Teff = $Tref $Teff_min $Teff_max


set model1 = `ls ${model_path}/p${Teff_min}_g${logg_min}*z${met_min}* | xargs -n1 basename | grep -v contopac`
 echo $model1
set model2 = `ls ${model_path}/p${Teff_min}_g${logg_min}*z${met_max}* | xargs -n1 basename | grep -v contopac`
 echo $model2
set model3 = `ls ${model_path}/p${Teff_min}_g${logg_max}*z${met_min}* | xargs -n1 basename | grep -v contopac`
 echo $model3
set model4 = `ls ${model_path}/p${Teff_min}_g${logg_max}*z${met_max}* | xargs -n1 basename | grep -v contopac`
 echo $model4
set model5 = `ls ${model_path}/p${Teff_max}_g${logg_min}*z${met_min}* | xargs -n1 basename | grep -v contopac`
 echo $model5
set model6 = `ls ${model_path}/p${Teff_max}_g${logg_min}*z${met_max}* | xargs -n1 basename | grep -v contopac`
 echo $model6
set model7 = `ls ${model_path}/p${Teff_max}_g${logg_max}*z${met_min}* | xargs -n1 basename | grep -v contopac`
 echo $model7
set model8 = `ls ${model_path}/p${Teff_max}_g${logg_max}*z${met_max}* | xargs -n1 basename | grep -v contopac`
 echo $model8

#enter here the values requested for the interpolated model

echo $model

if ($?model) then
  set modele_out = ${model}
else
  set modele_out =  ${Tref}g${loggref}m0.0z${zref}.int
endif

echo $modele_out


#format of output: 'TurboSpectrum',  'Atlas', 'MOOG',default
if ($?outtype) then
  echo $outtype
else
  set outtype = 'TurboSpectrum'
endif

#### the test option is set to .true. if you want to plot comparison model (model_test)
set test = '.false.'
set model_test = '/home/masseron/SPECTRUM/models/marcs_package/ftp/sun.mod'
#'MARCSbin', 'MARCSweb', 'Atlas' or 'Phoenix' format
set model_test_type = 'MARCSweb'
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
'${modele_out}'
${modtype}
${Tref}
${loggref}
${zref}
${outtype}
${test}
$model_test
$model_test_type
EOF
