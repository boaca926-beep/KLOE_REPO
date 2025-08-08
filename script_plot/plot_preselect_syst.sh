#!/bin/bash

cut_indx=3

## Initialization
CUT_LABEL=("egammamin"
	   "Zvmax"
	   "Rhovmax"
	   "nb_sigma_T_clust"
	  )

CUT_TITLE=("E^{min}_{clust} [MeV]" 
	   "Z^{max}_{v} [cm]"
	   "#rho^{max}_{v} [cm]"
	   "n [#sigma_{t}] [ps]"
	  )

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

header="../header/plot_omega_syst.h"
sed -i 's|\(const TString syst_path =\)\(.*\)|\1 "'"${syst_path}"'";|' "$header"
sed -i 's|\(const TString cut_label =\)\(.*\)|\1 "'"${cut_label}"'";|' "$header"
sed -i 's|\(const TString cut_title =\)\(.*\)|\1 "'"${cut_title}"'";|' "$header"

# Syst. main folder
plot_folder=../../plot_results/plot_preselect_syst_$cut_label/
if [[ -d "$plot_folder" ]]; then
    ls $plot_folder
    echo "Remove ${plot_folder}"
    rm -rf $plot_folder
else
    echo "Create ${plot_folder}"
fi
mkdir $plot_folder
ls ../../plot_results/

# Path
#plot_header=../header/plot_sfw_syst.h
#echo -e 'const TString plot_folder = "";' > $plot_header
#echo -e 'const TString outputPlot = "";' >> $plot_header


# Log files 
log=${plot_folder}log_plot.txt
touch $log

PARA_LABEL=("omega_mass" "omega_width" "BB") #list of parameter name
PARA_TITLE=("M_{#omega} [MeV\/c^{2}]" "#Gamma_{#omega} [MeV]" "B_{ee}B_{3#pi}") #list of parameter title
PARA_UNIT=("" "" "" "" "" "") #list of parameter unit

for ((i=0; i<${#PARA_LABEL[@]}; ++i)); do

    para_label=${PARA_LABEL[i]}
    para_title=${PARA_TITLE[i]}
    para_unit=${PARA_UNIT[i]}
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
    echo '  gROOT->ProcessLine(".L ../run_plot/plot_omega_syst.C");' >> $plot_script
    echo '  gROOT->ProcessLine("plot_omega_syst()");' >> $plot_script
    echo '}' >> $plot_script
    root -l -n -q -b $plot_script >> $log
    
done

rm $plot_script

merge_script=merge_script.C
echo '#include <iostream>' > $merge_script
echo "void merge_script() {" >> $merge_script
echo '  gROOT->ProcessLine(".L ../run_plot/mergePlots.C");' >> $merge_script
echo '  gROOT->ProcessLine("mergePlots()");' >> $merge_script
echo '}' >> $merge_script
root -l -n -q -b $merge_script

rm $merge_script

## Move result folder to
#dir=../../plot_results/plot_preselect_syst_${cut_label}

#if [[ -d "$dir" ]]; then
#    echo remove $dir
#    rm -rf $dir
#else
#    echo "${dir} does not exist!"
#    mkdir $dir
#    #exit 1
#fi

#cp -rf $plot_folder ../../plot_results/
#ls ../../plot_results/
