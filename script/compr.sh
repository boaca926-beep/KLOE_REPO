#!/bin/bash

# Main folder
compr_path=../../results_compr
if [[ -d "$compr_path" ]]; then
    ls $compr_path
    echo "Remove $compr_path"
    rm -rf $compr_path
else
    echo "Create ${compr_path}"
fi
mkdir ${compr_path}
echo "Results comparison folder is created at ${compr_path}"
ls ${compr_path}

path_header=../header/compr.h
sed -i 's|\(const TString outputFile =\)\(.*\)|\1 "'"${compr_path}"'";|' "$path_header"

mass_script=mass_script.C
echo '#include <iostream>' > $mass_script
echo "void mass_script() {" >> $mass_script
echo '  gROOT->ProcessLine(".L ../run/mass_compr.C");' >> $mass_script
echo '  gROOT->ProcessLine("mass_compr()");' >> $mass_script

#echo '  gROOT->ProcessLine(".L width_mass_new.C");' >> $mass_script
#echo '  gROOT->ProcessLine("width_compr()");' >> $mass_script

#echo '  gROOT->ProcessLine(".L ../run/branch_compr.C");' >> $mass_script
#echo '  gROOT->ProcessLine("branch_compr()");' >> $mass_script
echo '}' >> $mass_script
root -l -n -q -b $mass_script

branch_script=branch_script.C
echo '#include <iostream>' > $branch_script
echo "void branch_script() {" >> $branch_script
echo '  gROOT->ProcessLine(".L ../run/branch_compr.C");' >> $branch_script
echo '  gROOT->ProcessLine("branch_compr()");' >> $branch_script
echo '}' >> $branch_script
root -l -n -q -b $branch_script

rm $mass_script
rm $branch_script
