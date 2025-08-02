#!/bin/bash

# background rejecion 0: chi2; 1: angle; 2: deltaE; 3: beta
# fit 4: mass_sigma_nb; 6: fit_range
# normalzation 5: sfw2d_sigma_nb
# luminosity 7: lumi_nb

cut_indx=3
cut_type=select_cut 

## Initialization
CUT_LABEL=("chi2_cut" "angle_cut" "deltaE_cut" "beta_cut" "mass_sigma_nb" "sfw2d_sigma_nb" "fit_range" "lumi_nb")
CUT_TITLE=("#chi^{2}_{7C}" "#angle_{#gamma#gamma} [#circ]" "E_{miss} [MeV]" "c_{1}" "mass bin width variation [#deltaM_{3#pi}]" "#DeltaE_{#gamma_{3}}#times#DeltaM_{2#pi} bin width variation [%]" "Fitting range variation [#deltaM_{3#pi}]" "#DeltaL_{int} [#deltaL_{int}/L_{int}]") 
ERR_TYPE=(1 2 2 2 2 2 2 2) 

cut_label=${CUT_LABEL[$cut_indx]}
cut_title=${CUT_TITLE[$cut_indx]}
err_type=${ERR_TYPE[$cut_indx]}
echo "cut_label=${cut_label}; cut_title=${cut_title}; err_type=$err_type"

sample_type=norm # mini norm
#syst_path=../../result_syst/${sample_type}_${cut_label}/
syst_path=../../select_cut_syst_${sample_type}_TDATA_${cut_label}/

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
sed -i 's/\(int err_type =\)\(.*\)/\1 '$err_type';/' "$header"
    
# Syst. main folder
plot_folder=../../plot_results/plot_${sample_type}_${cut_label}/
if [[ -d "$plot_folder" ]]; then
    ls $plot_folder
    echo "Remove ${plot_folder}"
    rm -rf $plot_folder
else
    echo "Create ${plot_folder}"
fi
mkdir $plot_folder
ls ../../plot_results/

# Log files 
log=${plot_folder}log_plot.txt
touch $log

PARA_LABEL=("omega_mass" "omega_width" "BB") #list of parameter name
PARA_TITLE=("M_{#omega} [MeV\/c^{2}]" "#Gamma_{#omega} [MeV]" "B_{ee}B_{3#pi}") #list of parameter title

#PARA_LABEL=("omega_mass") #list of parameter name
#PARA_TITLE=("M_{#omega} [MeV\/c^{2}]") #list of parameter title

for ((j=0; j<${#PARA_LABEL[@]}; ++j)); do

    para_label=${PARA_LABEL[j]}
    para_title=${PARA_TITLE[j]}
    
    echo $j $para_label $para_title

    outputPlot=${plot_folder}
    
    sed -i 's|\(const TString outputPlot =\)\(.*\)|\1 "'"${outputPlot}"'";|' "$header"
    sed -i 's/\(const int para_indx =\)\(.*\)/\1 '$j';/' "$header"
    sed -i 's|\(const TString para_label =\)\(.*\)|\1 "'"${para_label}"'";|' "$header"
    sed -i 's|\(const TString para_title =\)\(.*\)|\1 "'"${para_title}"'";|' "$header"

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
#dir=../../plot_results/plot_${sample_type}_${cut_label}

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
