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

compr_script=compr_script.C
echo '#include <iostream>' > $compr_script
echo "void compr_script() {" >> $compr_script
echo '  gROOT->ProcessLine(".L ../run/mass_compr.C");' >> $compr_script
echo '  gROOT->ProcessLine("mass_compr()");' >> $compr_script

#echo '  gROOT->ProcessLine(".L width_compr_new.C");' >> $compr_script
#echo "  gROOT->ProcessLine(width_compr_new());" >> $compr_script

#echo '  gROOT->ProcessLine(".L branch_compr_new.C");' >> $compr_script
#echo "  gROOT->ProcessLine(branch_compr_new());" >> $compr_script
echo '}' >> $compr_script
root -l -n -q -b $compr_script
rm $compr_script
