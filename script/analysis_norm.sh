#!/bin/bash

## Specify DATA
sample_size=norm #norm; small; mini; chain
exp_type=TDATA 
gsf=1 

input_path=../../input_${sample_size}_${exp_type}
result_path=../../analysis_${sample_size}_${exp_type}

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
#echo "" > ${log_omega_fit}
touch $log_omega_fit

# initialize path.h
path_header=../header/path.h

echo -e 'const TString outputCut = "";' > $path_header
echo -e 'const TString outputHist = "";' >> $path_header
echo -e 'const TString outputOmega = "";' >> $path_header
echo -e 'const TString outputSfw2D = "";' >> $path_header
echo -e 'const TString outputSfw1D = "";' >> $path_header
echo -e 'const TString exp_type = "'$exp_type'";' >> $path_header
echo -e "double gsf = $gsf;" >> $path_header

outputCut=${input_path}/cut/
outputHist=${input_path}/hist/
outputSfw2D=${input_path}/sfw2d/
outputSfw1D=${input_path}/sfw1d/

sed -i 's|\(const TString outputCut =\)\(.*\)|\1 "'"${outputCut}"'";|' "$path_header"
sed -i 's|\(const TString outputHist =\)\(.*\)|\1 "'"${outputHist}"'";|' "$path_header"
sed -i 's|\(const TString outputOmega =\)\(.*\)|\1 "'"${omega_path}"'";|' "$path_header"

sed -i 's|\(const TString outputSfw2D =\)\(.*\)|\1 "'"${outputSfw2D}"'";|' "$path_header"
sed -i 's|\(const TString outputSfw1D =\)\(.*\)|\1 "'"${outputSfw1D}"'";|' "$path_header"

sed -i 's|\(const TString exp_type =\)\(.*\)|\1 "'"${exp_type}"'";|' "$path_header"
sed -i 's|\(double gsf =\)\(.*\)|\1 '$gsf';|' "$path_header"

# initialize sfw1d.txt and sfw2d.txt
#sfw1d_path=${input_path}/sfw1d/sfw1d.txt
#sfw2d_path=${input_path}/sfw2d/sfw2d.txt

#echo ${sfw2d_path} ${sfw1d_path}

#if [[ -f "${sfw1d_path}" ]]; then
    #echo "Remove ${sfw1d_path} ..."
    #rm -rf ${sfw1d_path}
    #exit 1
#else
    #echo "Folder ${sfw1d_path} doesn't exist!"
    #exit 1
#fi    

#
omega_fit_script=omega_fit_script.C
echo '#include <iostream>' > $omega_fit_script
echo "void omega_fit_script() {" >> $omega_fit_script
echo 'gROOT->ProcessLine(".L ../run/omega_fit.C");' >> $omega_fit_script
echo 'gROOT->ProcessLine("omega_fit()");' >> $omega_fit_script
echo '}' >> $omega_fit_script
root -l -n -q -b $omega_fit_script >> ${log_omega_fit}
echo "Omega parameters are extracted!"

rm $omega_fit_script
