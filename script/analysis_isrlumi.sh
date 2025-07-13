#!/bin/bash

## Specify DATA or UFO
sample_size=norm # norm; small; mini

exp_type=TDATA # DATA
gsf=1 # DATA

result_path=../../input_isrlumi_${sample_size}_${exp_type}

if [[ -d "$result_path" ]]; then
    #echo "Remove ${result_path} ..."
    rm -rf $result_path
fi    

omega_path=${result_path}/omega_fit/
log_path=${result_path}/log/

mkdir ${result_path} # result folder
mkdir ${omega_path} # omega parameters
mkdir ${log_path} # log files   

log_omega_fit=${log_path}log_omega_fit.txt
echo "" > ${log_omega_fit}
#touch $log_omega_fit

path_header=../header/path.h

outputCut=../../input_norm_TDATA/cut/
outputHist=../../input_norm_TDATA/hist/

sed -i 's|\(const TString outputCut =\)\(.*\)|\1 "'"${outputCut}"'";|' "$path_header"
sed -i 's|\(const TString outputHist =\)\(.*\)|\1 "'"${outputHist}"'";|' "$path_header"
sed -i 's|\(const TString outputOmega =\)\(.*\)|\1 "'"${omega_path}"'";|' "$path_header"
sed -i 's|\(const TString exp_type =\)\(.*\)|\1 "'"${exp_type}"'";|' "$path_header"
sed -i 's|\(double gsf =\)\(.*\)|\1 '$gsf';|' "$path_header"

omega_fit_script=omega_fit_script.C
echo '#include <iostream>' > $omega_fit_script
echo "void omega_fit_script() {" >> $omega_fit_script
echo 'gROOT->ProcessLine(".L ../run/omega_fit_isrlumi.C");' >> $omega_fit_script
echo 'gROOT->ProcessLine("omega_fit_isrlumi()");' >> $omega_fit_script
echo '}' >> $omega_fit_script
root -l -n -q -b $omega_fit_script >> ${log_omega_fit}
echo "Omega parameters are extracted!"

rm $omega_fit_script
