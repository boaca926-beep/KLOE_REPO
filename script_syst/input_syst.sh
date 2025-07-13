#!/bin/bash

sample_size=norm # norm; small; mini
vari_indx=3 # Variable index: 0: egammamin 1: Zvmax 2: Rhovmax 3: nb_sigma_T_clust
step_nb=1 # Number of variations (one side)

exp_type=TDATA # TDATA 
gsf=1 # TDATA

## Parameters
VARI=("egammamin" "Zvmax" "Rhovmax" "nb_sigma_T_clust") # Labels
NORM=(15 10 4 3) # Norminial value
STEP_SIZE=(5 0.1 0.1 0.5) # Variation step size

variable=${VARI[${vari_indx}]}
vari_norm=${NORM[${vari_indx}]} 
vari_step=${STEP_SIZE[${vari_indx}]}
vari_min=$(echo "$vari_norm-$step_nb*$vari_step" | bc)
vari_max=$(echo "$vari_norm+$step_nb*$vari_step" | bc)

echo "Variable: ${variable}"
echo "Range [$vari_min, $vari_max]"
echo "Variation step size: $vari_step" 
echo "Number of variations (one side): $step_nb"
echo "Sample size $sample_size"

## Initialize the normial conditions
egammamin=15
Rhovmax=4
Zvmax=10
nb_sigma_T_clust=3

class_header=../header/MyClass.h
sed -i 's/\(egammamin =\)\(.*\)/\1 '$egammamin';/' $class_header 
sed -i 's/\(Rhovmax =\)\(.*\)/\1 '$Rhovmax';/' $class_header
sed -i 's/\(Zvmax =\)\(.*\)/\1 '$Zvmax';/' $class_header   
sed -i 's/\(nb_sigma_T_clust =\)\(.*\)/\1 '$nb_sigma_T_clust';/' $class_header   

# tree_cut
chi2_cut=43
angle_cut=138
deltaE_cut=-150
beta_cut=1.98
c0=0.11
c1=0.8
cut_value=0

cut_header=../header/cut_para.h
echo -e 'const double chi2_cut = -1;' > $cut_header
echo -e 'const double angle_cut = -1;' >> $cut_header
echo -e 'const double deltaE_cut = -1;' >> $cut_header
echo -e 'const double beta_cut = -1;' >> $cut_header
echo -e 'const double c0 = -1;' >> $cut_header
echo -e 'const double c1 = -1;' >> $cut_header
echo -e 'const TString cut_nm = "";' >> $cut_header
echo -e 'double cut_value = -1;' >> $cut_header

sed -i 's/\(const double chi2_cut =\)\(.*\)/\1 '$chi2_cut';/' $cut_header
sed -i 's/\(const double angle_cut =\)\(.*\)/\1 '$angle_cut';/' $cut_header
sed -i 's/\(const double deltaE_cut =\)\(.*\)/\1 '$deltaE_cut';/' $cut_header
sed -i 's/\(const double beta_cut =\)\(.*\)/\1 '$beta_cut';/' $cut_header
sed -i 's/\(const double c0 =\)\(.*\)/\1 '$c0';/' $cut_header
sed -i 's/\(const double c1 =\)\(.*\)/\1 '$c1';/' $cut_header

sed -i 's/\(const TString cut_nm =\)\(.*\)/\1 "'$variable'";/' $cut_header
sed -i 's/\(double cut_value =\)\(.*\)/\1 '$cut_value';/' $cut_header

# histo
mass_sigma_nb=1
sfw2d_sigma_nb=1

hist_header=../header/hist.h
sed -i 's/\(const double mass_sigma_nb =\)\(.*\)/\1 '$mass_sigma_nb';/' $hist_header
sed -i 's/\(const double sfw2d_sigma_nb =\)\(.*\)/\1 '$sfw2d_sigma_nb';/' $hist_header

# omega_fit
fit_min=760
fit_max=800

omega_header=../header/omega_fit.h
sed -i 's/\(const double fit_min =\)\(.*\)/\1 '$fit_min';/' $omega_header
sed -i 's/\(const double fit_max =\)\(.*\)/\1 '$fit_max';/' $omega_header

# sm_para
Lumi_tot=1724470

sm_header=../header/sm_para.h
sed -i 's/\(const double Lumi_tot =\)\(.*\)/\1 '$Lumi_tot';/' $sm_header

# Syst. main folder
syst_path=../../presel_syst_$exp_type_$variable
if [[ -d "$syst_path" ]]; then
    ls $syst_path
    echo "Remove syst. folder"
    rm -rf $syst_path
else
    echo "Create ${syst_path}"
fi
mkdir ${syst_path}
echo "Syst. folder is created at ${syst_path}"
ls ${syst_path}

# Path
path_header=../header/path.h
echo -e 'const TString rootFile = "";' > $path_header
echo -e 'const TString sampleFile = "";' >> $path_header
echo -e 'const TString outputCut = "";' >> $path_header
echo -e 'const TString sig_path = "";' >> $path_header
echo -e 'const TString outputGen = "";' >> $path_header
echo -e 'const TString outputHist = "";' >> $path_header
echo -e 'const TString outputSfw2D = "";' >> $path_header
echo -e 'const TString outputSfw1D = "";' >> $path_header
echo -e 'const TString outputOmega = "";' >> $path_header
#echo -e 'const TString outputPath = "";' >> $path_header
echo -e 'const TString data_type = "";' >> $path_header
echo -e 'const TString exp_type = "'$exp_type'";' >> $path_header
echo -e "double gsf = $gsf;" >> $path_header

# Store output path
path_output=${syst_path}/path_output.txt
touch $path_output

sfw_output=${syst_path}/sfw_output.txt
touch $sfw_output

#sed -i 's|\(const TString outputPath =\)\(.*\)|\1 "'"${path_output}"'";|' "$path_header"

## Samples
DATA_TYPE=("sig" "ksl" "exp" "eeg" "ufo")
#sample_path=../../path_$sample_size # path_norm; path_small
sample_path=../path_${sample_size}/ 

## Start loop

for (( j=0; j<=step_nb*2; j ++))
    
do min=$(echo "$vari_min+$j*$vari_step" | bc)

   # Result folders
   result_path=$syst_path/result$j
   input_path=${result_path}/input/
   cut_path=${result_path}/cut/
   gen_path=${result_path}/gen/
   hist_path=${result_path}/hist/
   sfw2d_path=${result_path}/sfw2d/
   sfw1d_path=${result_path}/sfw1d/
   omega_path=${result_path}/omega_fit/
   log_path=${result_path}/log/
   
   if [[ -d "$result_path" ]]; then
       #echo "Remove ${result_path} ..."
       rm -rf $result_path
   fi    
   
   mkdir ${result_path} # result folder
   mkdir ${input_path} # input root files
   mkdir ${cut_path} # trees: after all cuts
   mkdir ${gen_path} # tree: signal MC generated
   mkdir ${hist_path} # histos
   mkdir ${sfw2d_path} # mc normalization
   mkdir ${sfw1d_path} # mc signal tuning
   mkdir ${omega_path} # omega parameters
   mkdir ${log_path} # log files   
   echo "Results folder is created at ${result_path}"

   # Log files 
   log_input=${log_path}log_input.txt
   #echo "" > ${log_input}
   touch $log_input
   
   log_cut=${log_path}log_cut.txt
   #echo "" > ${log_cut}
   touch $log_cut
   
   log_omega_fit=${log_path}log_omega_fit.txt
   #echo "" > ${log_omega_fit}
   touch $log_omega_fit
   echo "Initializing $path_header and log files!"

   log_hist=${log_path}log_hist.txt
   #echo "" > ${log_hist}
   touch $log_hist

   log_sfw2d=${log_path}log_sfw2d.txt
   #echo "" > ${log_sfw2d}
   touch $log_sfw2d
   
   log_sfw1d=${log_path}log_sfw1d.txt
   #echo "" > ${log_sfw1d}
   touch $log_sfw1d
   
   log_omega_fit=${log_path}log_omega_fit.txt
   #echo "" > ${log_omega_fit}
   touch $log_omega_fit

   ## Variations
   echo "Variation value: $min"
   sed -i 's/\('$variable' =\)\(.*\)/\1 '$min';/' ../header/MyClass.h
   sed -i 's/\(double cut_value =\)\(.*\)/\1 '$min';/' ../header/cut_para.h
   
   # loop over sample types
   for ((i=0;i<${#DATA_TYPE[@]};++i)); do

       ## create and update tree.h
       data_type=${DATA_TYPE[i]}
       
       INPUT_FILE=${sample_path}${data_type}_path
       ROOT_FILE=${result_path}/input/${data_type}
       echo $INPUT_FILE
       #echo $ROOT_FILE
       
       sed -i 's|\(const TString rootFile =\)\(.*\)|\1 "'"${INPUT_FILE}"'";|' "$path_header"
       sed -i 's|\(const TString sampleFile =\)\(.*\)|\1 "'"${ROOT_FILE}"'";|' "$path_header"
       
       ## create input trees
       class_script=class_script.C
       echo "void class_script() {" > $class_script
       echo '  gROOT->ProcessLine(".L ../run/MyClass.C");' >> $class_script
       echo '  gROOT->ProcessLine(".L ../run/Analys_class.C");' >> $class_script
       echo '  gROOT->ProcessLine("Analys_class(rootFile, sampleFile)");' >> $class_script
       echo '}' >> $class_script
       root -l -n -q -b $class_script >> ${log_input}
       
       ## Selection cuts
       sed -i 's|\(const TString outputCut =\)\(.*\)|\1 "'"${cut_path}"'";|' "$path_header"
       sed -i 's|\(const TString data_type =\)\(.*\)|\1 "'"${data_type}"'";|' "$path_header"
    
       tree_cut_script=tree_cut_script.C
       echo '#include <iostream>' > $tree_cut_script
       echo "void tree_cut_script() {" >> $tree_cut_script
       echo 'gROOT->ProcessLine(".L ../run/tree_cut.C");' >> $tree_cut_script
       echo 'gROOT->ProcessLine("tree_cut()");' >> $tree_cut_script
       echo '}' >> $tree_cut_script
       root -l -n -q -b $tree_cut_script >> ${log_cut}
   done
   echo "Selection cuts applied!"

   rm $class_script
   rm $tree_cut_script
   
   ## Signal MC generated
   sed -i 's|\(const TString sig_path =\)\(.*\)|\1 "'"${input_path}"'";|' "$path_header"
   sed -i 's|\(const TString outputGen =\)\(.*\)|\1 "'"${gen_path}"'";|' "$path_header"
   
   tree_gen_script=tree_gen_script.C
   echo '#include <iostream>' > $tree_gen_script
   echo "void tree_gen_script() {" >> $tree_gen_script
   echo 'gROOT->ProcessLine(".L ../run/tree_gen.C");' >> $tree_gen_script
   echo 'gROOT->ProcessLine("tree_gen()");' >> $tree_gen_script
   echo '}' >> $tree_gen_script
   root -l -n -q -b $tree_gen_script
   echo "Signal MC is generated!"

   ## Histos
   sed -i 's|\(const TString outputHist =\)\(.*\)|\1 "'"${hist_path}"'";|' "$path_header"
   hist_script=hist_script.C
   echo '#include <iostream>' > $hist_script
   echo "void hist_script() {" >> $hist_script
   echo 'gROOT->ProcessLine(".L ../run/gethist.C");' >> $hist_script
   echo 'gROOT->ProcessLine("gethist()");' >> $hist_script
   echo '}' >> $hist_script
   root -l -n -q -b $hist_script >> ${log_hist}
   echo "Histos are created!"

   ## Normalization
   sed -i 's|\(const TString outputSfw2D =\)\(.*\)|\1 "'"${sfw2d_path}"'";|' "$path_header"
   sfw2d_script=sfw2d_script.C
   echo '#include <iostream>' > $sfw2d_script
   echo "void sfw2d_script() {" >> $sfw2d_script
   echo 'gROOT->ProcessLine(".L ../run/sfw2d.C");' >> $sfw2d_script
   echo 'gROOT->ProcessLine("sfw2d()");' >> $sfw2d_script
   echo '}' >> $sfw2d_script
   root -l -n -q -b $sfw2d_script >> ${log_sfw2d}
   echo "MC normalization!"
   
   ## MC signal tuning
   sed -i 's|\(const TString outputSfw1D =\)\(.*\)|\1 "'"${sfw1d_path}"'";|' "$path_header"
   sfw1d_script=sfw1d_script.C
   echo '#include <iostream>' > $sfw1d_script
   echo "void sfw1d_script() {" >> $sfw1d_script
   echo 'gROOT->ProcessLine(".L ../run/sfw1d.C");' >> $sfw1d_script
   echo 'gROOT->ProcessLine("sfw1d()");' >> $sfw1d_script
   echo '}' >> $sfw1d_script
   root -l -n -q -b $sfw1d_script >> ${log_sfw1d}
   echo "MC signal tuning!"
   
   ## Omega parameters
   sed -i 's|\(const TString outputOmega =\)\(.*\)|\1 "'"${omega_path}"'";|' "$path_header"
   omega_fit_script=omega_fit_script.C
   echo '#include <iostream>' > $omega_fit_script
   echo "void omega_fit_script() {" >> $omega_fit_script
   echo 'gROOT->ProcessLine(".L ../run/omega_fit.C");' >> $omega_fit_script
   echo 'gROOT->ProcessLine("omega_fit()");' >> $omega_fit_script
   echo '}' >> $omega_fit_script
   root -l -n -q -b $omega_fit_script >> ${log_omega_fit}
   echo "Omega parameters are extracted!"
   
   ## Cut file
   #echo "${cut_path}tree_cut.root" >> $path_output # store result_path
   echo "${sfw2d_path}sfw2d.root" >> $sfw_output
   echo "${omega_path}omega_fit.root" >> $path_output
   
done

echo "Output path"
cat $path_output

## Remove script
rm $tree_gen_script
rm $hist_script
rm $sfw2d_script
rm $sfw1d_script
rm $omega_fit_script
   
   








