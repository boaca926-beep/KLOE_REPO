#!/bin/bash

variable="fit_range"
step_nb=4 # Number of variations (one side)
sample_type=norm #norm, chain, mini

## Initialize (nominal conditions)
# input
egammamin=15
Rhovmax=4
Zvmax=10
nb_sigma_T_clust=3

# tree_cut
chi2_cut=43
angle_cut=138
deltaE_cut=-150
beta_cut=1.98
c0=0.11
c1=0.8
cut_nm=""
cut_value=0

# histo
mass_sigma_nb=1
sfw2d_sigma_nb=1

# omega_fit
fit_min_sigma_nb=0;
fit_max_sigma_nb=0;

fit_min=760;
fit_max=800;

# sm_para
Lumi_tot=1724470

# update header files
class_header=../header/MyClass.h
sed -i 's/\(egammamin =\)\(.*\)/\1 '$egammamin';/' $class_header
sed -i 's/\(Rhovmax =\)\(.*\)/\1 '$Rhovmax';/' $class_header
sed -i 's/\(Zvmax =\)\(.*\)/\1 '$Zvmax';/' $class_header   
sed -i 's/\(nb_sigma_T_clust =\)\(.*\)/\1 '$nb_sigma_T_clust';/' $class_header   

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

hist_header=../header/hist.h
sed -i 's/\(const double mass_sigma_nb =\)\(.*\)/\1 '$mass_sigma_nb';/' $hist_header
sed -i 's/\(const double sfw2d_sigma_nb =\)\(.*\)/\1 '$sfw2d_sigma_nb';/' $hist_header

omega_header=../header/omega_fit.h
sed -i 's/\(const double fit_min =\)\(.*\)/\1 '$fit_min';/' $omega_header
sed -i 's/\(const double fit_max =\)\(.*\)/\1 '$fit_max';/' $omega_header

sm_header=../header/sm_para.h
sed -i 's/\(const double Lumi_tot =\)\(.*\)/\1 '$Lumi_tot';/' $sm_header

path_header=../header/path.h
exp_type=TDATA 
echo -e 'const TString sampleFile = "";' > $path_header
#echo -e 'const TString sig_path = "";' >> $path_header
echo -e 'const TString outputGen = "";' >> $path_header
echo -e 'const TString outputCut = "";' >> $path_header
echo -e 'const TString outputHist = "";' >> $path_header
echo -e 'const TString outputSfw2D = "";' >> $path_header
echo -e 'const TString outputSfw1D = "";' >> $path_header
echo -e 'const TString outputOmega = "";' >> $path_header
echo -e 'const TString data_type = "";' >> $path_header
echo -e 'const TString exp_type = "'$exp_type'";' >> $path_header
echo -e "double gsf = 1;" >> $path_header

## Variation parameters
vari_step=0.5
vari_norm=0 

VARI1=("fit_min_sigma_nb") # Labels
VARI2=("fit_max_sigma_nb")

vari1_min=$(echo "$vari_norm-$step_nb*$vari_step" | bc)
vari1_max=$(echo "$vari_norm+$step_nb*$vari_step" | bc)

vari2_min=$(echo "$vari_norm-$step_nb*$vari_step" | bc)
vari2_max=$(echo "$vari_norm+$step_nb*$vari_step" | bc)

echo "Variable1: ${VARI1[0]}, range [$vari1_min, $vari1_max], step size: $vari_step, # variations (one side): $step_nb"
echo "Variable2: ${VARI2[0]}, range [$vari2_min, $vari2_max], step size: $vari_step, # variations (one side): $step_nb"

## Check syst. input folder
norm_path=../../result_input/input_${sample_type}_TDATA
#norm_path=../../result/input_${sample_type}_TDATA
input_path=${norm_path}/input/

if [[ -d "$norm_path" ]]; then
    echo input_path: $input_path
    ls $norm_path/input
else
    echo "${norm_path} does not exist!"
    exit 1
fi

#sed -i 's|\(const TString sig_path =\)\(.*\)|\1 "'"${input_path}"'";|' "$path_header"
sed -i 's|\(const TString outputGen =\)\(.*\)|\1 "'"${norm_path}/gen/"'";|' "$path_header"

## Syst. main folder
syst_path=../../result_syst/${sample_type}_$variable

#: <<'END_COMMENT'
if [[ -d "$syst_path" ]]; then
    echo "Remove syst. folder"
    rm -rf $syst_path
else
    echo "Create ${syst_path}"
fi
#END_COMMENT

mkdir ${syst_path}
ls ${syst_path}
echo "Syst. folder is created at ${syst_path}"

## Log and txt files
log_path=${syst_path}/log/
mkdir ${log_path} # log files   

#log_gen=${log_path}gen
#mkdir $log_gen

log_cut=${log_path}cut
mkdir $log_cut

log_hist=${log_path}hist
mkdir $log_hist

log_sfw2d=${log_path}sfw2d
mkdir $log_sfw2d

log_sfw1d=${log_path}sfw1d
mkdir $log_sfw1d

log_omega=${log_path}omega_fit
mkdir $log_omega

# 
path_output=${syst_path}/path_output.txt
touch $path_output

sfw2d_path_output=${syst_path}/sfw2d_path_output.txt
touch $sfw2d_path_output

## Create sub-folders
cut_path=${syst_path}/cut/
mkdir ${cut_path} # tree: after all cuts

#gen_path=${syst_path}/gen/
#mkdir ${gen_path} # tree: signal MC generated

hist_path=${syst_path}/hist/
mkdir ${hist_path} # histos

sfw2d_path=${syst_path}/sfw2d/
mkdir ${sfw2d_path} # mc normalization

sfw1d_path=${syst_path}/sfw1d/
mkdir ${sfw1d_path} # mc signal tuning

omega_path=${syst_path}/omega_fit/
mkdir ${omega_path} # omega parameters

## Signal MC generated
#log_gen_file=$log_gen/log.txt
#touch $log_gen_file
   
#tree_gen_script=tree_gen_script.C
#echo '#include <iostream>' > $tree_gen_script
#echo "void tree_gen_script() {" >> $tree_gen_script
#echo 'gROOT->ProcessLine(".L ../run/tree_gen.C");' >> $tree_gen_script
#echo 'gROOT->ProcessLine("tree_gen()");' >> $tree_gen_script
#echo '}' >> $tree_gen_script
#root -l -n -q -b $tree_gen_script >> ${log_gen_file}
#echo "Signal MC is generated!"

## Start loop
echo "fit range nominal=[$fit_min $fit_max]"

cut_header=../header/cut_para.h

for (( j=0; j<=step_nb*2; j ++))
    
do min=$(echo "${fit_min}+($vari1_min+$j*$vari_step)*2.65" | bc)
   max=$(echo "${fit_max}-($vari2_min+$j*$vari_step)*2.65" | bc)

   ## Variations
   #sed -i 's/\(const double fit_min =\)\(.*\)/\1 '$fit_min';/' $omega_header
   #sed -i 's/\(const double fit_max =\)\(.*\)/\1 '$fit_max';/' $omega_header

   indx=$(echo "$j+1" | bc)

   sed -i 's/\(const double fit_min =\)\(.*\)/\1 '$min';/' $omega_header
   sed -i 's/\(const double fit_max =\)\(.*\)/\1 '$max';/' $omega_header
   
   echo "index=$indx; ${VARI1[0]}, ${VARI2[0]}=$(echo "$vari1_min+$j*$vari_step" | bc); fit range=[$min, $max]"
   sed -i 's/\(const TString cut_nm =\)\(.*\)/\1 "'$variable'";/' $cut_header
   sed -i 's/\(double cut_value =\)\(.*\)/\1 '$(echo "$vari1_min+$j*$vari_step" | bc)';/' $cut_header
   
   tree_subfolder=${cut_path}tree$j/ 
   mkdir $tree_subfolder # sub folder in cut/tree_cut0, 1, 2 ...
   echo $tree_subfolder

   hist_subfolder=${hist_path}hist$j/ 
   mkdir $hist_subfolder # sub folder in hist/hist0, 1, 2 ...
   echo $hist_subfolder
   
   sed -i 's|\(const TString outputCut =\)\(.*\)|\1 "'"${tree_subfolder}"'";|' "$path_header"

   sfw2d_subfolder=${sfw2d_path}sfw2d$j/ 
   mkdir $sfw2d_subfolder # sub folder in sfw2d/sfw2d0, 1, 2 ...
   echo $sfw2d_subfolder
   echo "${sfw2d_subfolder}sfw2d.root" >> $sfw2d_path_output 
   
   sfw1d_subfolder=${sfw1d_path}sfw1d$j/ 
   mkdir $sfw1d_subfolder # sub folder in sfw1d/sfw1d0, 1, 2 ...
   echo $sfw1d_subfolder

   omega_subfolder=${omega_path}omega_fit$j/ 
   mkdir $omega_subfolder # sub folder in omega/omega_fit0, 1, 2 ...
   echo $omega_subfolder
   echo "${omega_subfolder}omega_fit.root" >> $path_output 
   
   # log files
   log_cut_file=$log_cut/log$j.txt
   touch $log_cut_file

   log_hist_file=$log_hist/log$j.txt
   touch $log_hist_file

   log_sfw2d_file=$log_sfw2d/log$j.txt
   touch $log_sfw2d_file

   log_sfw1d_file=$log_sfw1d/log$j.txt
   touch $log_sfw1d_file

   log_omega_file=$log_omega/log$j.txt
   touch $log_omega_file

   DATA_TYPE=("sig" "ksl" "exp" "eeg" "ufo")

   # loop over sample types
   for ((i=0;i<${#DATA_TYPE[@]};++i)); do

       ## Create trees after selection cuts
       data_type=${DATA_TYPE[i]}
       
       ROOT_FILE=${input_path}${data_type}
       echo $ROOT_FILE

       sed -i 's|\(const TString sampleFile =\)\(.*\)|\1 "'"${ROOT_FILE}"'";|' "$path_header"
       sed -i 's|\(const TString data_type =\)\(.*\)|\1 "'"${data_type}"'";|' "$path_header"

       tree_cut_script=tree_cut_script.C
       echo '#include <iostream>' > $tree_cut_script
       echo "void tree_cut_script() {" >> $tree_cut_script
       echo 'gROOT->ProcessLine(".L ../run/tree_cut.C");' >> $tree_cut_script
       echo 'gROOT->ProcessLine("tree_cut()");' >> $tree_cut_script
       echo '}' >> $tree_cut_script
       root -l -n -q -b $tree_cut_script >> ${log_cut_file}

       rm ${tree_cut_script}

   done
   echo "Selection cuts applied!"

   ## Histos
   sed -i 's|\(const TString outputHist =\)\(.*\)|\1 "'"${hist_subfolder}"'";|' "$path_header"
   hist_script=hist_script.C
   echo '#include <iostream>' > $hist_script
   echo "void hist_script() {" >> $hist_script
   echo 'gROOT->ProcessLine(".L ../run/gethist.C");' >> $hist_script
   echo 'gROOT->ProcessLine("gethist()");' >> $hist_script
   echo '}' >> $hist_script
   root -l -n -q -b $hist_script >> ${log_hist_file}
   echo "Histos are created!"

   ## Normalization
   sed -i 's|\(const TString outputSfw2D =\)\(.*\)|\1 "'"${sfw2d_subfolder}"'";|' "$path_header"
   sfw2d_script=sfw2d_script.C
   echo '#include <iostream>' > $sfw2d_script
   echo "void sfw2d_script() {" >> $sfw2d_script
   echo 'gROOT->ProcessLine(".L ../run/sfw2d.C");' >> $sfw2d_script
   echo 'gROOT->ProcessLine("sfw2d()");' >> $sfw2d_script
   echo '}' >> $sfw2d_script
   root -l -n -q -b $sfw2d_script >> ${log_sfw2d_file}
   echo "MC normalization!"

   ## MC signal tuning
   sed -i 's|\(const TString outputSfw1D =\)\(.*\)|\1 "'"${sfw1d_subfolder}"'";|' "$path_header"
   sfw1d_script=sfw1d_script.C
   echo '#include <iostream>' > $sfw1d_script
   echo "void sfw1d_script() {" >> $sfw1d_script
   echo 'gROOT->ProcessLine(".L ../run/sfw1d.C");' >> $sfw1d_script
   echo 'gROOT->ProcessLine("sfw1d()");' >> $sfw1d_script
   echo '}' >> $sfw1d_script
   root -l -n -q -b $sfw1d_script >> ${log_sfw1d_file}
   echo "MC signal tuning!"

   ## Omega parameters
   sed -i 's|\(const TString outputOmega =\)\(.*\)|\1 "'"${omega_subfolder}"'";|' "$path_header"
   omega_fit_script=omega_fit_script.C
   echo '#include <iostream>' > $omega_fit_script
   echo "void omega_fit_script() {" >> $omega_fit_script
   echo 'gROOT->ProcessLine(".L ../run/omega_fit.C");' >> $omega_fit_script
   echo 'gROOT->ProcessLine("omega_fit()");' >> $omega_fit_script
   echo '}' >> $omega_fit_script
   root -l -n -q -b $omega_fit_script >> ${log_omega_file}
   echo "Omega parameters are extracted!"

done

#rm ${tree_gen_script}
rm ${hist_script}
rm ${sfw2d_script}
rm ${sfw1d_script}
rm ${omega_fit_script}



