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

set model_path = models/MARCS_st_sph_t02_mod
set model_path = models/c+0.00/MARCS_st_sph_t02_mod

echo $model_path
#'MARCSbin', 'MARCSweb', 'Atlas' or 'Phoenix' format
set modtype = 'MARCSweb'

if (!(-s ${model_path}/grid.list)) then
  ls ${model_path}/s*g*m*z* | xargs -n1 basename > ${model_path}/grid.list
  echo grid.list created
endif

set model1 = s7250_g+2.0_m1.0_t02_x3_z-0.25_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod

 echo $model1
set model2 = s7250_g+2.0_m1.0_t02_x3_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod

 echo $model2
set model3 = s7250_g+2.5_m1.0_t02_x3_z-0.25_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod

 echo $model3
set model4 = s7250_g+2.5_m1.0_t02_x3_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod

 echo $model4
set model5 = s7500_g+2.0_m1.0_t02_x3_z-0.25_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod

 echo $model5
set model6 = s7500_g+2.0_m1.0_t02_x3_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod

 echo $model6
set model7 = s7500_g+2.5_m1.0_t02_x3_z-0.25_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod

 echo $model7
set model8 = s7500_g+2.5_m1.0_t02_x3_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod
 echo $model8

#enter here the values requested for the interpolated model 
if ($?model) then
  set modele_out = ${model}
else
  set modele_out =  ${Tref}g${loggref}m1.0z${zref}.int
endif

echo $modele_out


#format of output: 'TurboSpectrum',  'Atlas', 'MOOG',default
set outtype = 'TurboSpectrum'
#### the test option is set to .true. if you want to plot comparison model (model_test)
set test = '.true.'
set model_test = ${model_path}/s7000_g+1.5_m1.0_t02_x3_z-0.25_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod
echo $model_test
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
'$model_test'
$model_test_type
EOF

echo 'control plot loading...'
#---------------------------------------------------------------------------------------------------
# control plot with sm
if ( ${test} == ".true.") then
  set testsm = '1'
else
  set testsm = '0'
endif


plots/sm/sm -f /net/nas4/nebula/tmasseron/BACCHUS/plots/sm/.sm << eof 
data "models.sm" 
verbose 0
lines 3 0
read  {tau1 3 T1 4 Pe1 5 Pg1 6 taur1 2}
define ndp (dimen(tau1)-1)
define end (\$ndp+3)
do i=2,9 {
define beg (\$end+3)
define end (\$beg + \$ndp)
lines \$beg \$end
read  tau\$i 3 
read T\$i 4 
read Pe\$i 5 
read Pg\$i 6 
read taur\$i 2
}

if (${testsm}) {
define beg (\$end+3)
lines \$beg 0
read  {tau10 3 T10 4 Pe10 5 Pg10 6 taur10 2}
}

device postportfile interpol_check.ps
ctype green
toplabel T_{eff}=${Tref} logg=${loggref} z=${zref}
ctype black
location 3000 15000 17000 32000
ylabel log(Pe)
limits tau9 Pe9
box 0 2 0 0

lweight 4
ltype 0
ctype blue
relocate (4000 31000)
draw (4600 31000) 
putlabel 6 low T_{eff}
ctype red
relocate (4000 30500)
draw (4600 30500)
putlabel 6 up T_{eff}
ctype black
lweight 1
relocate (4000 30000)
draw (4600 30000)
putlabel 6 low logg
lweight 4
relocate (4000 29500)
draw (4600 29500)
putlabel 6 up logg
ltype 0
relocate (4000 29000)
draw (4600 29000)
putlabel 6 low [Fe/H]
ltype 1
relocate (4000 28500)
draw (4600 28500)
putlabel 6 up [Fe/H]

ctype blue
lweight 1
ltype 0
connect tau1 Pe1
ltype 1
connect tau2 Pe2
lweight 2
ltype 0
connect tau3 Pe3
ltype 1
connect tau4 Pe4
ctype red
lweight 1
ltype 0
connect tau5 Pe5
ltype 1
connect tau6 Pe6
lweight 2
ltype 0
connect tau7 Pe7
ltype 1
connect tau8 Pe8
ltype 0
lweight 1
ctype green
connect tau9 Pe9
ctype black
if (${testsm}) {
connect tau10 Pe10   
}

location 3000 15000 3000 17000
xlabel log (\tau_{5000})
ylabel T
limits tau9 T9
box
ctype blue
connect tau1 T1
ltype 1
connect tau2 T2
lweight 2
ltype 0
connect tau3 T3
ltype 1
connect tau4 T4
ctype red
lweight 1
ltype 0
connect tau5 T5
ltype 1
connect tau6 T6
ltype 0
lweight 2
connect tau7 T7
ltype 1
connect tau8 T8
ltype 0
ctype green
connect tau9 T9
ctype black
if (${testsm}) {
connect tau10 T10  
}
lweight 1

location 19000 32000 17000 32000
xlabel log (\tau_{5000})
ylabel log(Pg)
limits tau9 Pg9
box  
ctype blue
connect tau1 Pg1
ltype 1
connect tau2 Pg2
ltype 0
lweight 2
connect tau3 Pg3
ltype 1
connect tau4 Pg4
ctype red
ltype 0
lweight 1
connect tau5 Pg5
ltype 1
connect tau6 Pg6
ltype 0
lweight 2
connect tau7 Pg7
ltype 1
connect tau8 Pg8
ltype 0
ctype green
connect tau9 Pg9
ctype black
if (${testsm}) {
connect tau10 Pg10 
}
lweight 1

if (${testsm}) {
location 19300 32000 3000 11000
ylabel actual error (%)
location 20000 32000 3000 6000
spline tau9  T9 tau10 Tbis
limits tau10 (abs(Tbis-T10)/T10 * 100)
connect  tau10 (abs(Tbis-T10)/T10 * 100)
box
ylabel T
xlabel log (\tau_{5000})

location 20000 32000 6000 9000
spline tau9 Pe9 tau10 Pebis
limits tau10 (abs(10**Pebis-10**Pe10)/10**Pe10 * 100)
box 0 2 0 0
connect  tau10  (abs(10**Pebis-10**Pe10)/10**Pe10 * 100)
ylabel Pe

location 20000 32000 9000 12000
spline tau9 Pg9 tau10 Pgbis
limits tau10 (abs(10**Pgbis-10**Pg10)/10**Pg10 * 100)
connect  tau10 (abs(10**Pgbis-10**Pg10)/10**Pg10 * 100)
box 0 2 0 0
ylabel Pg
}

eof


#gv -orientation=portrait interpol_check.ps&

#rm -f models.sm

