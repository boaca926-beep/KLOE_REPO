#!/bin/bash

## Test .C files: sfw2d.C, sfw1d.C and omega_fit.C
## Check consistency with ../../crx3pi_UFOKslCut/output/

## Test sfw2d.C

## Output folders
#run_test=../../run_test
#sfw2d=${run_test}/sfw2d/
#sfw1d=${run_test}/sfw1d/
#omega_fit=${run_test}/omega_fit/
#log=${run_test}/log/

#if [[ -d "$run_test" ]]; then
#    #echo "Remove ${run_test} ..."
#    rm -rf $run_test
#fi

#mkdir ${run_test} # result folder
#mkdir ${sfw2d} # sfw2d folder
#mkdir ${omega_fit} # omega folder
#mkdir ${log} # log folder

cut_path="../../crx3pi_UFOKslCut/output/"
hist_path="../../crx3pi_UFOKslCut/output/"
code_ref="../code_ref/"

# Path
path_header=../header/path.h
echo -e 'const TString outputCut = "";' > $path_header
echo -e 'const TString outputHist = "";' >> $path_header
echo -e 'const TString outputSfw2D = "";' >> $path_header
echo -e 'const TString outputSfw1D = "";' >> $path_header
echo -e 'const TString outputOmega = "";' >> $path_header

sed -i 's|\(const TString outputCut =\)\(.*\)|\1 "'"${cut_path}"'";|' "$path_header"
sed -i 's|\(const TString outputHist =\)\(.*\)|\1 "'"${hist_path}"'";|' "$path_header"
sed -i 's|\(const TString outputSfw2D =\)\(.*\)|\1 "'"${code_ref}"'";|' "$path_header"
sed -i 's|\(const TString outputSfw1D =\)\(.*\)|\1 "'"${code_ref}"'";|' "$path_header"
sed -i 's|\(const TString outputOmega =\)\(.*\)|\1 "'"${code_ref}"'";|' "$path_header"

# Log files 
log_sfw2d=${code_ref}log_sfw2d.txt
#echo "" > ${log_sfw2d}
touch $log_sfw2d

log_sfw1d=${code_ref}log_sfw1d.txt
#echo "" > ${log_sfw1d}
touch $log_sfw1d

log_omega=${code_ref}log_omega.txt
#echo "" > ${log_omega}
touch $log_omega

## Normalization
sfw2d_script=sfw2d_script.C
echo '#include <iostream>' > $sfw2d_script
echo "void sfw2d_script() {" >> $sfw2d_script
echo 'gROOT->ProcessLine(".L ../run/sfw2d.C");' >> $sfw2d_script
echo 'gROOT->ProcessLine("sfw2d()");' >> $sfw2d_script
echo '}' >> $sfw2d_script
#root -l -n -q -b $sfw2d_script >> ${log_sfw2d}
echo "MC normalization!"

sfw1d_script=sfw1d_script.C
echo '#include <iostream>' > $sfw1d_script
echo "void sfw1d_script() {" >> $sfw1d_script
echo 'gROOT->ProcessLine(".L ../run/sfw1d.C");' >> $sfw1d_script
echo 'gROOT->ProcessLine("sfw1d()");' >> $sfw1d_script
echo '}' >> $sfw1d_script
#root -l -n -q -b $sfw1d_script >> ${log_sfw1d}
echo "MC signal tuning!"

omega_fit_script=omega_fit_script.C
echo '#include <iostream>' > $omega_fit_script
echo "void omega_fit_script() {" >> $omega_fit_script
echo 'gROOT->ProcessLine(".L ../run/omega_fit.C");' >> $omega_fit_script
echo 'gROOT->ProcessLine("omega_fit()");' >> $omega_fit_script
echo '}' >> $omega_fit_script
root -l -n -q -b $omega_fit_script >> ${log_omega}
echo "Omega parameters are extracted!"

# Remove script
rm $sfw2d_script
rm $sfw1d_script
rm $omega_fit_script

