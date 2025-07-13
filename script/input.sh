#!/bin/bash

## define output folder: smaple size, type of DATA ana folder name
sample_size=chain #chain norm
exp_type=TDATA
gsf=1 # DATA

result_path=../../input_${sample_size}_${exp_type}
input_path=${result_path}/input/
sample_path=../path_${sample_size}/
cut_path=${result_path}/cut/
log_path=${result_path}/log/

## Create main folder
#: <<'END_COMMENT'
# create the main folder
if [[ -d "$result_path" ]]; then
    #echo "Remove ${result_path} ..."
    rm -rf $result_path
fi    
#END_COMMENT

mkdir ${result_path} # result folder
mkdir ${input_path} # input root files
mkdir ${cut_path} # trees: after all cuts
mkdir ${log_path} # log files   

## Initialize the normial conditions
egammamin=15
Rhovmax=4
Zvmax=10
nb_sigma_T_clust=3
chi2_cut=43
angle_cut=138
deltaE_cut=-150
beta_cut=1.98
cut_value=3.8
c0=0.11
c1=0.8
   
sed -i 's/\(egammamin =\)\(.*\)/\1 '$egammamin';/' ../header/MyClass.h
sed -i 's/\(Rhovmax =\)\(.*\)/\1 '$Rhovmax';/' ../header/MyClass.h
sed -i 's/\(Zvmax =\)\(.*\)/\1 '$Zvmax';/' ../header/MyClass.h   

cut_header=../header/cut_para.h
sed -i 's/\(const double chi2_cut =\)\(.*\)/\1 '$chi2_cut';/' $cut_header
sed -i 's/\(const double angle_cut =\)\(.*\)/\1 '$angle_cut';/' $cut_header
sed -i 's/\(const double deltaE_cut =\)\(.*\)/\1 '$deltaE_cut';/' $cut_header
sed -i 's/\(const double beta_cut =\)\(.*\)/\1 '$beta_cut';/' $cut_header
sed -i 's/\(const double c0 =\)\(.*\)/\1 '$c0';/' $cut_header
sed -i 's/\(const double c1 =\)\(.*\)/\1 '$c1';/' $cut_header


## Log files 
log_input=${log_path}log_input.txt
#echo "" > ${log_input}
touch $log_input

log_cut=${log_path}log_cut.txt
#echo "" > ${log_cut}
touch $log_cut

## create sub folder store output files
DATA_TYPE=("sig" "ksl" "exp" "eeg" "ufo")
#DATA_TYPE=("sig")

path_header=../header/path.h
echo -e 'const TString rootFile = "";' > $path_header
echo -e 'const TString sampleFile = "";' >> $path_header
echo -e 'const TString inputFile = "";' >> $path_header
echo -e 'const TString chainFile = "";' >> $path_header
echo -e 'const TString outputCut = "";' >> $path_header
echo -e 'const TString sig_path = "";' >> $path_header
echo -e 'const TString outputGen = "";' >> $path_header
echo -e 'const TString outputHist = "";' >> $path_header
echo -e 'const TString outputSfw2D = "";' >> $path_header
echo -e 'const TString outputSfw1D = "";' >> $path_header
echo -e 'const TString outputOmega = "";' >> $path_header
echo -e 'const TString data_type = "";' >> $path_header
echo -e 'const TString exp_type = "'$exp_type'";' >> $path_header
echo -e "double gsf = $gsf;" >> $path_header

for ((i=0;i<${#DATA_TYPE[@]};++i)); do

    data_type=${DATA_TYPE[i]}
    #echo ${data_type}
    
    if [[ -d "${result_path}/${DATA_TYPE[i]}" ]]; then
	echo "Remove ${DATA_TYPE[i]} ..."
	rm -rf ${result_path}/${DATA_TYPE[i]}
    fi    

    ROOT_PATH=${input_path}${data_type}

    mkdir ${ROOT_PATH}
    
    sed -i 's|\(const TString data_type =\)\(.*\)|\1 "'"${data_type}"'";|' "$path_header"
    
    # loop over path 1, 2 and 3
    for ((k=1;k<=3;++k)); do

	INPUT_FILE=${sample_path}${data_type}_path$k
	echo $INPUT_FILE

	ROOT_FILE=${ROOT_PATH}/${data_type}${k}
	echo $ROOT_FILE
    
	sed -i 's|\(const TString rootFile =\)\(.*\)|\1 "'"${INPUT_FILE}"'";|' "$path_header"
	sed -i 's|\(const TString sampleFile =\)\(.*\)|\1 "'"${ROOT_FILE}"'";|' "$path_header"
	
	## Input trees
	run_script=run_script.C
	
	echo "void run_script() {" > $run_script
	echo '  gROOT->ProcessLine(".L ../run/MyClass.C");' >> $run_script
	echo '  gROOT->ProcessLine(".L ../run/Analys_class.C");' >> $run_script
	echo '  gROOT->ProcessLine("Analys_class(rootFile, sampleFile)");' >> $run_script
	echo '}' >> $run_script
	root -l -n -q -b $run_script >> ${log_input}
 	
    done

    ## Merge input root files
    sed -i 's|\(const TString inputFile =\)\(.*\)|\1 "'"${input_path}${data_type}"'";|' "$path_header"
    sed -i 's|\(const TString chainFile =\)\(.*\)|\1 "'"${input_path}"'";|' "$path_header"

    chains_script=chains_script.C
	
    echo "void chains_script() {" > $chains_script
    echo '  gROOT->ProcessLine(".L ../run/getchains.C");' >> $chains_script
    echo '  gROOT->ProcessLine("getchains(inputFile, chainFile)");' >> $chains_script
    echo '}' >> $chains_script
    root -l -n -q -b $chains_script
    
    ls ${ROOT_PATH}
    rm -rf ${ROOT_PATH}

    echo "${data_type} input merged!"

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

## Signal MC generated
gen_path=${result_path}/gen/
mkdir ${gen_path} # tree: signal MC generated

log_gen=${log_path}log_gen.txt
touch $log_gen

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
hist_path=${result_path}/hist/
mkdir ${hist_path} # histos

log_hist=${log_path}log_hist.txt
touch $log_hist

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
sfw2d_path=${result_path}/sfw2d/
sfw1d_path=${result_path}/sfw1d/

mkdir ${sfw2d_path} # mc normalization
mkdir ${sfw1d_path} # mc signal tuning

log_sfw2d=${log_path}log_sfw2d.txt
touch $log_sfw2d

log_sfw1d=${log_path}log_sfw1d.txt
touch $log_sfw1d

sed -i 's|\(const TString outputSfw2D =\)\(.*\)|\1 "'"${sfw2d_path}"'";|' "$path_header"
sfw2d_script=sfw2d_script.C
sed -i 's|\(const TString outputSfw1D =\)\(.*\)|\1 "'"${sfw1d_path}"'";|' "$path_header"
sfw1d_script=sfw1d_script.C

echo '#include <iostream>' > $sfw2d_script
echo "void sfw2d_script() {" >> $sfw2d_script
echo 'gROOT->ProcessLine(".L ../run/sfw2d.C");' >> $sfw2d_script
echo 'gROOT->ProcessLine("sfw2d()");' >> $sfw2d_script
echo '}' >> $sfw2d_script
root -l -n -q -b $sfw2d_script >> ${log_sfw2d}
echo "MC normalization!"

## MC signal tuning
echo '#include <iostream>' > $sfw1d_script
echo "void sfw1d_script() {" >> $sfw1d_script
echo 'gROOT->ProcessLine(".L ../run/sfw1d.C");' >> $sfw1d_script
echo 'gROOT->ProcessLine("sfw1d()");' >> $sfw1d_script
echo '}' >> $sfw1d_script
root -l -n -q -b $sfw1d_script >> ${log_sfw1d}
echo "MC signal tuning!"

## Omega parameters
omega_path=${result_path}/omega_fit/
mkdir ${omega_path} # omega parameters
log_omega_fit=${log_path}log_omega_fit.txt
touch $log_omega_fit
sed -i 's|\(const TString outputOmega =\)\(.*\)|\1 "'"${omega_path}"'";|' "$path_header"

omega_fit_script=omega_fit_script.C
echo '#include <iostream>' > $omega_fit_script
echo "void omega_fit_script() {" >> $omega_fit_script
echo 'gROOT->ProcessLine(".L ../run/omega_fit.C");' >> $omega_fit_script
echo 'gROOT->ProcessLine("omega_fit()");' >> $omega_fit_script
echo '}' >> $omega_fit_script
root -l -n -q -b $omega_fit_script >> ${log_omega_fit}
echo "Omega parameters are extracted!"

## Remove script
rm ${run_script}
rm ${tree_cut_script}
rm ${tree_gen_script}
rm ${hist_script}
rm ${sfw2d_script}
rm ${sfw1d_script}
rm ${omega_fit_script}
