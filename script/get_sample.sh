#!/bin/bash

# Samples
#DATA_TYPE=("sig" "ksl" "exp" "eeg" "ufo")
DATA_TYPE=("sig" "ksl" "exp" "eeg")

cut_path=../../
input_path=../../kloe_sample
#input_path=../../input_norm_TDATA

if [[ -d "$input_path" ]]; then
    echo "${input_path} exists"
    #rm -rf $input_path
fi    

cut_file="${cut_path}tree_pre.root"
if [[ -f "$cut_file" ]]; then
    echo "$cut_file exists"
    rm -rf $cut_file
else
    echo "Error: $cut_file not found!"
    exit 1
fi    

mkdir ${input_path} # result folder

path_header=../header/path_sample.h
echo -e 'const TString data_type = "";' > $path_header
echo -e 'const TString sampleFile = "";' >> $path_header
echo -e 'const TString outputCut = "";' >> $path_header

echo "Looping over data samples ..."

input_path=${input_path}/input/
cut_path=${cut_path}

sed -i 's|\(const TString outputCut =\)\(.*\)|\1 "'"${cut_path}"'";|' "$path_header"
    
# Loop over data smaples
for ((i=0;i<${#DATA_TYPE[@]};++i)); do

    data_type=${DATA_TYPE[i]}
    
    INPUT_FILE=${input_path}${data_type}
    echo ${data_type} $INPUT_FILE

    sed -i 's|\(const TString data_type =\)\(.*\)|\1 "'"${data_type}"'";|' "$path_header"
    sed -i 's|\(const TString sampleFile =\)\(.*\)|\1 "'"${INPUT_FILE}"'";|' "$path_header"

    tree_sample_script=tree_sample_script.C
    echo '#include <iostream>' > $tree_sample_script
    echo "void tree_sample_script() {" >> $tree_sample_script
    echo 'gROOT->ProcessLine(".L ../run_vertex_bkgrej/tree_sample.C");' >> $tree_sample_script
    echo 'gROOT->ProcessLine("tree_sample()");' >> $tree_sample_script
    echo '}' >> $tree_sample_script
    root -l -n -q -b $tree_sample_script
    
done

rm $tree_sample_script

