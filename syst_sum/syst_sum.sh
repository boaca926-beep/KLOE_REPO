#!/bin/bash

#main_folder="../../SYST"
log_syst_summary=../../SYST/log_syst_summary.txt

if [[ -e "$log_syst_summary" ]]; then
    echo update $log_syst_summary
    rm -rf $log_syst_summary
else
    echo "" > $log_syst_summary
    #touch $log_syst_summary
    #exit 1
fi

syst_summary_script=syst_summary_script.C
echo '#include <iostream>' > $syst_summary_script
echo "void syst_summary_script() {" >> $syst_summary_script
echo 'gROOT->ProcessLine(".L ../run/syst_summary.C");' >> $syst_summary_script
echo 'gROOT->ProcessLine("syst_summary()");' >> $syst_summary_script
echo '}' >> $syst_summary_script
root -l -n -q -b $syst_summary_script >> ${log_syst_summary}
echo "Total systematic uncertainty is obtained!"

rm $syst_summary_script
