#!/bin/bash

## run > ./hist.sh ANALYSIS2025

## specify data types: sig ksl exp eeg ufo
DATA_TYPE=("sig" "ksl" "exp" "eeg" "ufo")
#DATA_TYPE=("sig")

input_path=/media/bo/8E97-E8DD/${1}

cut_header=cut.h
    

## create and update analysis folder
if [[ -d "$input_path/cut" ]]; then
    echo update ${input_path} ...
    rm ${input_path}/cut/*.root
    rm ${input_path}/hist/*.root
    rm ${input_path}/mc_norm/*.root

else
    mkdir ${input_path}/cut
    mkdir ${input_path}/hist
    mkdir ${input_path}/mc_norm
fi

log_cut_file=${input_path}/log_cut.txt
echo "" > ${log_cut_file}


for ((i=0;i<${#DATA_TYPE[@]};++i)); do

    ## create and update tree.h
    data_type=${DATA_TYPE[i]}
    #echo ${data_type}

    sed -i 's|\(const TString data_type =\)\(.*\)|\1 "'"${data_type}"'";|' "$cut_header"
    
    sed -i 's|\(const TString cut_in =\)\(.*\)|\1 "'"${input_path}/input/${data_type}.root"'";|' "$cut_header"
    sed -i 's|\(const TString cut_out =\)\(.*\)|\1 "'"${input_path}/cut/tree_cut.root"'";|' "$cut_header"
    
    ## cut trees
    tree_cut_script=tree_cut_script.C
    echo '#include <iostream>' > $tree_cut_script
    echo "void tree_cut_script() {" >> $tree_cut_script
    echo 'gROOT->ProcessLine(".L tree_cut.C");' >> $tree_cut_script
    echo 'gROOT->ProcessLine("tree_cut()");' >> $tree_cut_script
    echo '}' >> $tree_cut_script
    root -l -n -q -b $tree_cut_script >> ${log_cut_file}
    
done

hist_header=hist.h
sed -i 's|\(const TString hist_out =\)\(.*\)|\1 "'"${input_path}/hist/hist.root"'";|' "$hist_header"

echo "Preparing hist.root ..."
hist_script=hist_script.C
echo '#include <iostream>' > $hist_script
echo "void hist_script() {" >> $hist_script
echo 'gROOT->ProcessLine(".L gethist.C");' >> $hist_script
echo 'gROOT->ProcessLine("gethist()");' >> $hist_script
echo '}' >> $hist_script
root -l -n -q -b $hist_script

## Remove script files
rm $tree_cut_script
rm $hist_script
