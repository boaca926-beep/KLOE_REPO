#!/bin/bash

input_path=/media/bo/8E97-E8DD/ANALYSIS2025

## create and update mc normalization folder
if [[ -d "$input_path/omega_para" ]]; then
    echo update ${input_path} ...
    rm ${input_path}/omega_para/*.root
else
    mkdir ${input_path}/omega_para
fi

log_analys_file=${input_path}/log_analys.txt
echo "" > ${log_analys_file}

echo "Obtain omeaga parameters ..."
omega_para_header=omega_para.h

sed -i 's|\(const TString hist_in =\)\(.*\)|\1 "'"${input_path}/hist/hist.root"'";|' "$omega_para_header"
sed -i 's|\(const TString cut_in =\)\(.*\)|\1 "'"${input_path}/cut"'";|' "$omega_para_header"
sed -i 's|\(const TString omega_para_out =\)\(.*\)|\1 "'"${input_path}/omega_para"'";|' "$omega_para_header"


log_omega_para_file=${input_path}/log_omega_para.txt
echo "" > ${log_omega_para_file}

omega_para_script=omega_para_script.C
echo '#include <iostream>' > $omega_para_script
echo "void omega_para_script() {" >> $omega_para_script
echo 'gROOT->ProcessLine(".L omega_para.C");' >> $omega_para_script
echo 'gROOT->ProcessLine("omega_para()");' >> $omega_para_script
echo '}' >> $omega_para_script
root -l -n -q -b $omega_para_script >> ${log_omega_para_file}

rm $omega_para_script
