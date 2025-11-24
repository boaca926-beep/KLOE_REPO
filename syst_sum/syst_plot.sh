#!/bin/bash

para_indx=2
PARA_NAME=("Mass_omega" "Gam_omega" "BB") #list of parameter name
#PARA_NAME=("Mass_omega")
#PARA_NAME=("Gam_omega")
#PARA_NAME=("BB")

input=../header/syst_plot.h

for ((i=0; i<${#PARA_NAME[@]}; ++i)); do


    #echo ${PARA_NAME[i]}

    sed -i 's/\(const int para_indx =\)\(.*\)/\1 '$para_indx';/' $input
    sed -i 's/\(const TString para_name =\)\(.*\)/\1 "'${PARA_NAME[$para_indx]}'";/' $input

    syst_script=syst_script.C
    echo '#include <iostream>' > $syst_script
    echo "void syst_script() {" >> $syst_script
    echo '  gROOT->ProcessLine(".L syst_plot.C");' >> $syst_script
    echo '  gROOT->ProcessLine("syst_plot()");' >> $syst_script
    echo '}' >> $syst_script
    #root -l -n -q -b $syst_script 
    
done

rm $syst_script
