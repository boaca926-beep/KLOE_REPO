#!/bin/bash

echo "Preparing tree_gen.root ..."
run_script2=run_script2.C
echo '#include <iostream>' > $run_script2
echo "void run_script2() {" >> $run_script2
echo 'gROOT->ProcessLine(".L tree_gen.C");' >> $run_script2
echo 'gROOT->ProcessLine("tree_gen()");' >> $run_script2
echo '}' >> $run_script2
#root -l -n -q -b $run_script2
rm $run_script2

## Histos
hist_script=hist_script.C
echo '#include <iostream>' > $hist_script
echo "void hist_script() {" >> $hist_script
echo 'gROOT->ProcessLine(".L gethist.C");' >> $hist_script
echo 'gROOT->ProcessLine("gethist()");' >> $hist_script
echo '}' >> $hist_script
#root -l -n -q -b $hist_script 
echo "Histos are created!"

## Normalization
#sfw2d_script=sfw2d_script.C
#echo '#include <iostream>' > $sfw2d_script
#echo "void sfw2d_script() {" >> $sfw2d_script
#echo 'gROOT->ProcessLine(".L sfw2d.C");' >> $sfw2d_script
#echo 'gROOT->ProcessLine("sfw2d()");' >> $sfw2d_script
#echo '}' >> $sfw2d_script
#root -l -n -q -b $sfw2d_script
#echo "MC normalization!"

## MC signal tuning
#sfw1d_script=sfw1d_script.C
#echo '#include <iostream>' > $sfw1d_script
#echo "void sfw1d_script() {" >> $sfw1d_script
#echo 'gROOT->ProcessLine(".L sfw1d.C");' >> $sfw1d_script
#echo 'gROOT->ProcessLine("sfw1d()");' >> $sfw1d_script
#echo '}' >> $sfw1d_script
#root -l -n -q -b $sfw1d_script
#echo "MC signal tuning!"

log_omega_fit=log_omega_fit.txt
echo "" > ${log_omega_fit}
#touch $log_omega_fit 

## Omega parameters
omega_fit_script=omega_fit_script.C
echo '#include <iostream>' > $omega_fit_script
echo "void omega_fit_script() {" >> $omega_fit_script
echo 'gROOT->ProcessLine(".L omega_fit.C");' >> $omega_fit_script
echo 'gROOT->ProcessLine("omega_fit()");' >> $omega_fit_script
echo '}' >> $omega_fit_script
root -l -n -q -b $omega_fit_script >> $log_omega_fit
echo "Omega parameters are extracted!"

log_result=log_result.txt
touch $log_result 

result_script=result_script.C
echo '#include <iostream>' > $result_script
echo "void result_script() {" >> $result_script
echo 'gROOT->ProcessLine(".L result.C");' >> $result_script
echo 'gROOT->ProcessLine("result()");' >> $result_script
echo '}' >> $result_script
#root -l -n -q -b $result_script >> $log_result
echo "Omega parameters are extracted!"

rm $hist_script
rm $sfw2d_script
rm $sfw1d_script
rm $omega_fit_script
rm $result_script




