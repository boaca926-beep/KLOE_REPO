#!/bin/bash

PARA_NAME=("Mass_omega" "Gam_omega" "BB") #list of parameter name
PARA_INDEX=(0 1 2)

#PARA_NAME=("Mass_omega") #list of parameter name
#PARA_INDEX=(0)

input=../header/plot_models_syst.h

for ((i=0; i<${#PARA_NAME[@]}; ++i)); do

    #echo ${PARA_NAME[i]}

    sed -i 's/\(const TString para_name =\)\(.*\)/\1 "'${PARA_NAME[i]}'";/' $input
    sed -i 's/\(const int para_index =\)\(.*\)/\1 '${PARA_INDEX[i]}';/' $input
    

    syst_script=syst_script.C
    echo '#include <iostream>' > $syst_script
    echo "void syst_script() {" >> $syst_script
    echo '  gROOT->ProcessLine(".L ../run_plot/plot_models_syst.C");' >> $syst_script
    echo '  gROOT->ProcessLine("plot_models_syst()");' >> $syst_script
    echo '}' >> $syst_script
    root -l -n -q -b $syst_script 
    
done

rm $syst_script
