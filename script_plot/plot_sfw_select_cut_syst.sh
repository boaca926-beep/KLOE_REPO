#!/bin/bash

## Parameters
CUT_LABEL=("chi2_cut"
	   "angle_cut"
	   "deltaE_cut"
	  )

CUT_TITLE=("#chi^{2}_{7C}"
	   "#angle_{#gamma#gamma} [#circ]"
	   "E_{miss} [MeV]"
	  ) 

## Plot main folder
plot_folder=../../plot_sfw_select_cut_syst/

if [[ -d "$plot_folder" ]]; then
    echo "Remove the plot folder"
    rm -rf $plot_folder
else
    echo "Create ${plot_folder}"
fi
mkdir ${plot_folder}

header="../header/plot_sfw_syst.h"
sed -i 's|\(const TString pathType =\)\(.*\)|\1 "sfw2d_path_output.txt";|' "$header"
    

## Scaling factors
PARA_LABEL=("eeg_sfw" "isr3pi_sfw" "omegapi_sfw" "etagam_sfw" "ksl_sfw" "mcrest_sfw") #list of parameter name
PARA_TITLE=("ee#gamma" "signal" "#omega#pi" "#eta#gamma" "K_{S}K_{L}" "MC others") #list of parameter title
PARA_UNIT=("" "" "" "" "" "") #list of parameter unit

# loop over cuts 0: chi2 1: angle 2: deltaE 3: beta
for ((k=0;k<${#CUT_LABEL[@]};++k)); do

    cut_label=${CUT_LABEL[k]}
    cut_title=${CUT_TITLE[k]}
    echo ${cut_label} ${cut_title}


    echo $cut_label
    mkdir ${plot_folder}$cut_label

    outputPlot=${plot_folder}$cut_label/
  
    # check input files
    syst_path=../../select_cut_syst_norm_TDATA_${cut_label}/
    if [[ -d "$syst_path" ]]; then
	ls ${syst_path}sfw2d
    else
	echo "${syst_path} does not exist!"
	exit 1
    fi

    sed -i 's|\(const TString syst_path =\)\(.*\)|\1 "'"${syst_path}"'";|' "$header"
    sed -i 's|\(const TString output_path =\)\(.*\)|\1 "'"${plot_folder}"'";|' "$header"
    sed -i 's|\(const TString outputPlot =\)\(.*\)|\1 "'"${outputPlot}"'";|' "$header"

    sed -i 's|\(const TString cut_label =\)\(.*\)|\1 "'"${cut_label}"'";|' "$header"
    sed -i 's|\(const TString cut_title =\)\(.*\)|\1 "'"${cut_title}"'";|' "$header"
    
    # loop over sfw
    for ((i=0; i<${#PARA_LABEL[@]}; ++i)); do

	para_label=${PARA_LABEL[i]}
	para_title=${PARA_TITLE[i]}
	para_unit=${PARA_UNIT[i]}
	para_indx=$(echo "$i  " | bc)

	echo $para_indx $para_label $para_title $para_unit

	
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
	root -l -n -q -b $plot_script 
    	
    done
    
    #if (i==0); then
#	break
 #   fi
    
done

rm $plot_script

