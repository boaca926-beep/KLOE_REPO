#!/bin/bash

echo -e "\nPlotting histo comparison ..."

##################################################################
#VAR_NM=("deltaE")
#VAR_SYMB=("E_{diff}")
#UNIT=("[MeV]")                                                                                                                                                       
#XMIN=(-600) #-700, -500
#XMAX=(-50) #-200, 50
#BINS=(150) #150, 550

##################################################################
#VAR_NM=("IM_pi0_7C")
#VAR_SYMB=("M_{#gamma#gamma}")
#UNIT=("[MeV\/c^{2}]")                                                                    

#XMIN=(110)
#XMAX=(200)
#BINS=(180)

##################################################################
#VAR_NM=("angle_pi0gam12")
#VAR_SYMB=("#angle_{#gamma#gamma}")
#UNIT=("[#circ]")                                                                         

#XMIN=(0) #20
#XMAX=(180) #140
#BINS=(180) #120

##################################################################                        
#VAR_NM="betapi0"
#VAR_SYMB="#beta_{#pi}"
#UNIT=""                                                                                  

#XMIN=0.6
#XMAX=1
#BINS=150

##########################################################
#VAR_NM="lagvalue_min_7C"
#VAR_SYMB="#chi^{2}_{7C}"
#UNIT=""

#XMIN=0
#XMAX=50
#BINS=150

##########################################################
VAR_NM=("Eprompt_max")
VAR_SYMB=("E_{#gamma}^{max}")
UNIT=("[MeV]")

XMIN=(150)
XMAX=(350)
BINS=(150)

##########################################################
#VAR_NM=("Ephi_miss")
#VAR_SYMB=("E_{#phi-3#pi#gamma}^{miss}")
#UNIT=("[MeV]")

#XMIN=(-300)
#XMAX=(200)
#BINS=(150)

##########################################################
#VAR_NM=("IM3pi_7C")
#VAR_SYMB=("M_{3#pi}")
#UNIT=("[MeV]")

#XMIN=(520)
#XMAX=(580)
#BINS=(150)

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
    root -l -n -q -b $compr_script >> output.txt
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
