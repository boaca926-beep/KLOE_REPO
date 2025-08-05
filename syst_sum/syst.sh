#!/bin/bash

# preselection (all contri.)
# bkg rejection (all contr.: chi2_cut; angle_cut; deltaE_cut; beta_cut);
# fit (all contr.: mass_sigma_nb; fit_range);
# normalization (all contr.: sfw2d_sigma_nb)
# integrated lumi (lumi_nb)
#


CUT_TYPE=("egammamin" "nb_sigma_T_clust")
sample_type="preselect_syst"
syst_type="presel"

#CUT_TYPE=("chi2_cut" "angle_cut" "deltaE_cut" "beta_cut")
#sample_type="norm"
#syst_type="bkg"

#CUT_TYPE=("sfw2d_sigma_nb")
#sample_type="norm"
#syst_type="sfw2d"

#CUT_TYPE=("mass_sigma_nb" "fit_range")
#sample_type="norm"
#syst_type="fit"

#CUT_TYPE=("lumi_nb")
#sample_type="norm"
#syst_type="lumi"


main_folder="../../SYST/syst_${syst_type}"

if [[ -d "$main_folder" ]]; then
    echo main_folder: $main_folder
    rm -rf $main_folder
else
    echo "Create ${main_folder}!"
fi

mkdir ${main_folder}

PARA_LABEL=("omega_mass" "omega_width" "BB") #list of parameter name
PARA_TITLE=("M_{#omega} [MeV\/c^{2}]" "#Gamma_{#omega} [MeV]" "B_{ee}B_{3#pi}") #list of parameter title

for ((l=0;l<${#PARA_LABEL[@]};++l)); do

    para_label=${PARA_LABEL[l]}
    para_title=${PARA_TITLE[l]}
    
    echo $para_label $para_title

    input_folder=${main_folder}/${para_label}
    mkdir ${input_folder}
    echo "Syst. Summary folder is created at ${input_folder}"

    log_syst_summary=${input_folder}/log.txt
    touch $log_syst_summary

    input_path=${input_folder}/input.txt
    touch $input_path
    
    # update header files
    header=../header/syst_summary.h
    echo -e 'const TString input_folder = "";' > $header
    echo -e 'const TString para_label = "";' >> $header
    #echo -e 'const TString cut_type = "";' >> $header

    sed -i 's|\(const TString input_folder =\)\(.*\)|\1 "'"${input_folder}"'";|' "$header"
    sed -i 's|\(const TString para_label =\)\(.*\)|\1 "'"${para_label}"'";|' "$header"

    # loop cut
    for ((i=0;i<${#CUT_TYPE[@]};++i)); do

	cut_type=${CUT_TYPE[i]}
	DIR=../../plot_results/plot_${sample_type}_${cut_type}
	echo $DIR/${para_label}_${cut_type}.root
	#echo $cut_type
    
	# Loop through all .root files in the directory
	for file in "$DIR"/${para_label}_${cut_type}.root; do
	    
	    # Check if the file exists (handles case where no .root files are found)
	    if [ -e "$file" ]; then
		echo "Processing file: $file"
		echo "$file" >> "$input_path"
	    else
		exit 1
       	    fi
	
	done

    done
    # end the loop

    syst_summary_script=syst_summary_script.C
    echo '#include <iostream>' > $syst_summary_script
    echo "void syst_summary_script() {" >> $syst_summary_script
    echo 'gROOT->ProcessLine(".L ../run/syst.C");' >> $syst_summary_script
    echo 'gROOT->ProcessLine("syst()");' >> $syst_summary_script
    echo '}' >> $syst_summary_script
    root -l -n -q -b $syst_summary_script >> ${log_syst_summary}
    echo "Syst. are summerized!"
    
done

rm $syst_summary_script

