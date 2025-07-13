#!/bin/bash

log_syst_summary=../header/log_syst_summary.txt
echo "" > $log_syst_summary
#touch $log_syst_summary

syst_summary_script=syst_summary_script.C
echo '#include <iostream>' > $syst_summary_script
echo "void syst_summary_script() {" >> $syst_summary_script
echo 'gROOT->ProcessLine(".L ../run/syst_summary.C");' >> $syst_summary_script
echo 'gROOT->ProcessLine("syst_summary()");' >> $syst_summary_script
echo '}' >> $syst_summary_script
root -l -n -q -b $syst_summary_script >> ${log_syst_summary}
echo "Total systematic uncertainty is obtained!"

rm $syst_summary_script
