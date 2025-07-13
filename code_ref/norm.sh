#!/bin/bash

input_path=/media/bo/8E97-E8DD/ANALYSIS2025

## create and update mc normalization folder
if [[ -d "$input_path/mc_norm" ]]; then
    echo update ${input_path} ...
    rm ${input_path}/mc_norm/*.root
else
    mkdir ${input_path}/mc_norm
fi

log_norm_file=${input_path}/log_norm.txt
echo "" > ${log_norm_file}

echo "Preparing MC normalization and signal MC tuning ..."
sfw1d_header=sfw1d.h
sfw2d_header=sfw2d.h

sed -i 's|\(const TString hist_in =\)\(.*\)|\1 "'"${input_path}/hist/hist.root"'";|' "$sfw2d_header"
sed -i 's|\(const TString norm_out =\)\(.*\)|\1 "'"${input_path}/mc_norm"'";|' "$sfw2d_header"

sed -i 's|\(const TString hist_in =\)\(.*\)|\1 "'"${input_path}/hist/hist.root"'";|' "$sfw1d_header"
sed -i 's|\(const TString norm_out =\)\(.*\)|\1 "'"${input_path}/mc_norm"'";|' "$sfw1d_header"
sed -i 's|\(const TString cut_in =\)\(.*\)|\1 "'"${input_path}/cut"'";|' "$sfw1d_header"


sfw2d_script=sfw2d_script.C
echo '#include <iostream>' > $sfw2d_script
echo "void sfw2d_script() {" >> $sfw2d_script
echo 'gROOT->ProcessLine(".L sfw2d.C");' >> $sfw2d_script
echo 'gROOT->ProcessLine("sfw2d()");' >> $sfw2d_script
echo '}' >> $sfw2d_script
root -l -n -q -b $sfw2d_script >> ${log_norm_file}

sfw1d_script=sfw1d_script.C
echo '#include <iostream>' > $sfw1d_script
echo "void sfw1d_script() {" >> $sfw1d_script
echo 'gROOT->ProcessLine(".L sfw1d.C");' >> $sfw1d_script
echo 'gROOT->ProcessLine("sfw1d()");' >> $sfw1d_script
echo '}' >> $sfw1d_script
root -l -n -q -b $sfw1d_script >> ${log_norm_file}

rm $sfw1d_script
rm $sfw2d_script
