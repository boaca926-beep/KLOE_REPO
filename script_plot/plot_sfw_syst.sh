#!/bin/bash

cut_indx=0

## Initialization
CUT_LABEL=("egammamin"
	   "Zvmax"
	   "Rhovmax")

CUT_TITLE=("E^{min}_{clust} [MeV]" 
	   "Z^{max}_{v} [cm]"
	   "#rho^{max}_{v} [cm]")

cut_label=${CUT_LABEL[$cut_indx]}
cut_title=${CUT_TITLE[$cut_indx]}
echo ${cut_label} ${cut_title}

syst_path=../../presel_syst_$cut_label/
if [[ -d "$syst_path" ]]; then
    ls $syst_path
else
    echo "${syst_path} does not exist!"
    exit 1
fi

path_type="sfw_output.txt"
header="../header/plot_sfw_syst.h"
sed -i 's|\(const TString syst_path =\)\(.*\)|\1 "'"${syst_path}"'";|' "$header"
sed -i 's|\(const TString pathType =\)\(.*\)|\1 "'"${path_type}"'";|' "$header"
sed -i 's|\(const TString cut_label =\)\(.*\)|\1 "'"${cut_label}"'";|' "$header"
sed -i 's|\(const TString cut_title =\)\(.*\)|\1 "'"${cut_title}"'";|' "$header"

# Syst. main folder
plot_folder=../../plot_sfw_syst_$cut_label/
if [[ -d "$plot_folder" ]]; then
    ls $plot_folder
    echo "Remove plot folder"
    rm -rf $plot_folder
else
    echo "Create ${plot_folder}"
fi
mkdir $plot_folder

# Path
#plot_header=../header/plot_sfw_syst.h
#echo -e 'const TString plot_folder = "";' > $plot_header
#echo -e 'const TString outputPlot = "";' >> $plot_header


# Log files 
log=${plot_folder}log_plot.txt
touch $log


#para_indx=0
SFW_LABEL=("eeg_sfw" "isr3pi_sfw" "omegapi_sfw" "etagam_sfw" "ksl_sfw" "mcrest_sfw") #list of parameter name
SFW_TITLE=("ee#gamma" "signal" "#omega#pi" "#eta#gamma" "K_{S}K_{L}" "MC others") #list of parameter title
SFW_UNIT=("" "" "" "" "" "") #list of parameter unit

for ((i=0; i<${#SFW_LABEL[@]}; ++i)); do

    para_label=${SFW_LABEL[i]}
    para_title=${SFW_TITLE[i]}
    para_unit=${SFW_UNIT[i]}
    para_indx=$(echo "$i  " | bc)
    echo $para_indx $para_label $para_title $para_unit

    outputPlot=${plot_folder}
    #mkdir $outputPlot
    #echo "Plot folder is created at ${plot_folder}"

    sed -i 's|\(const TString outputPlot =\)\(.*\)|\1 "'"${outputPlot}"'";|' "$header"

    sed -i 's/\(const int para_indx =\)\(.*\)/\1 '$para_indx';/' "$header"
    sed -i 's|\(const TString para_label =\)\(.*\)|\1 "'"${para_label}"'";|' "$header"
    sed -i 's|\(const TString para_title =\)\(.*\)|\1 "'"${para_title} scaling factor"'";|' "$header"
    sed -i 's|\(const TString para_unit =\)\(.*\)|\1 "'"${para_unit}"'";|' "$header"
    

    plot_script=plot_script.C
    echo '#include <iostream>' > $plot_script
    echo "void plot_script() {" >> $plot_script
    echo '  gROOT->ProcessLine(".L ../run_plot/plot_sfw_syst.C");' >> $plot_script
    echo '  gROOT->ProcessLine("plot_sfw_syst()");' >> $plot_script
    echo '}' >> $plot_script
    root -l -n -q -b $plot_script >> $log
    
done

rm $plot_script
