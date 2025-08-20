#!/bin/bash

echo -e "\nPlotting histo comparison ..."

#VAR_NM="lagvalue_min_7C"                                                                                                                                            
#VAR_SYMB="#chi^{2}_{7C}"                                                                                                                                            
#UNIT=""

#XMIN=0                                                                                                                                                              
#XMAX=100                                                                                                                                                            
#BINS=200

##########################################################
VAR_NM=("Eprompt_max")
VAR_SYMB=("E_{#gamma}^{max}")
UNIT=("[MeV]")

XMIN=(150)
XMAX=(450)
BINS=(300)

output_folder="../../bkg_compr_"${VAR_NM[0]}

#check output folder and update output files
if [[ -d $output_folder ]]; then
    
    echo updating $output_folder
    #rm $output_folder/*.pdf
    #rm $output_folder/*.root
    
else
    
    echo root file $output_folder does not exsit;
    mkdir $output_folder
    
fi

path_header=../header/bkg_compr.h
sed -i 's|\(const TString output_folder =\)\(.*\)|\1 "'"${output_folder}"'";|' "$path_header"
    
for ((i=0;i<${#VAR_NM[@]};++i)); do

    sed -i 's/\(const int binsize =\)\(.*\)/\1 '${BINS[i]}';/' $path_header
    sed -i 's/\(const double var_min =\)\(.*\)/\1 '${XMIN[i]}';/' $path_header
    sed -i 's/\(const double var_max =\)\(.*\)/\1 '${XMAX[i]}';/' $path_header
    
    sed -i 's/\(const TString var_nm =\)\(.*\)/\1 "'${VAR_NM[i]}'";/' $path_header
    sed -i 's/\(const TString unit =\)\(.*\)/\1 "'${UNIT[i]}'";/' $path_header
    sed -i 's/\(const TString var_symb =\)\(.*\)/\1 "'${VAR_SYMB[i]}'";/' $path_header
    
    compr_script=compr_script.C
    echo '#include <iostream>' > $compr_script
    echo "void compr_script() {" >> $compr_script
    echo '  gROOT->ProcessLine(".L ../run/bkg_compr.C");' >> $compr_script
    echo '  gROOT->ProcessLine("bkg_compr()");' >> $compr_script
    echo '}' >> $compr_script
    #root -l -n -q -b $compr_script >> output.txt
    rm $compr_script

    plot_script=plot_script.C
    echo '#include <iostream>' > $plot_script
    echo "void plot_script() {" >> $plot_script
    echo '  gROOT->ProcessLine(".L ../run_plot/plot_compr.C");' >> $plot_script
    echo '  gROOT->ProcessLine("plot_compr()");' >> $plot_script
    echo '}' >> $plot_script
    root -l -n -q -b $plot_script >> output.txt
    rm $plot_script
    
    plot_script=plot_script.C
    echo '#include <iostream>' > $plot_script
    echo "void plot_script() {" >> $plot_script
    echo '  gROOT->ProcessLine(".L ../run_plot/plot_hist.C");' >> $plot_script
    echo '  gROOT->ProcessLine("plot_hist()");' >> $plot_script
    echo '}' >> $plot_script
    root -l -n -q -b $plot_script >> output.txt
    rm $plot_script
    
done
