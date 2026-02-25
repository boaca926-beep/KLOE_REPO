echo -e "\nPlotting histo comparison ..."

##################################################################
VAR_NM=("pull_E1" "pull_x1")
VAR_SYMB=("pull E1" "pull x1")
UNIT=("[MeV]" "")                                                                      
XMIN=(-5 -5) 
XMAX=(5 5) 
BINS=(200 200) 

output_folder="../../pulls_compr"

#check output folder and update output files
if [[ -d $output_folder ]]; then
    
    echo updating $output_folder
    #rm $output_folder/*.pdf
    rm $output_folder/*.root
    
else
    
    echo root file $output_folder does not exsit;
    mkdir $output_folder
    
fi

path_header=../header/pulls_compr.h
sed -i 's#\(const TString output_folder =\).*#\1 "'"${output_folder}"'";#' "$path_header"
#sed -i 's|\(const TString output_folder =\)\(.*\)|\1 "'"${output_folder}"'";|' "$path_header"
    
for ((i=0;i<${#VAR_NM[@]};++i)); do

    sed -i 's/\(const int binsize =\)\(.*\)/\1 '${BINS[i]}';/' $path_header
    sed -i 's/\(const double var_min =\)\(.*\)/\1 '${XMIN[i]}';/' $path_header
    sed -i 's/\(const double var_max =\)\(.*\)/\1 '${XMAX[i]}';/' $path_header
    
    sed -i 's/\(const TString var_nm =\)\(.*\)/\1 "'${VAR_NM[i]}'";/' $path_header
    sed -i 's/\(const TString unit =\)\(.*\)/\1 "'${UNIT[i]}'";/' $path_header
    #sed -i 's/\(const TString var_symb =\)\(.*\)/\1 "'${VAR_SYMB[i]}'";/' $path_header
    sed -i "s/\(const TString var_symb =\)\(.*\)/\1 \"${VAR_SYMB[i]}\";/" "$path_header"
    
    compr_script=compr_script.C
    echo '#include <iostream>' > $compr_script
    echo "void compr_script() {" >> $compr_script
    echo '  gROOT->ProcessLine(".L ../run/pulls_compr.C");' >> $compr_script
    echo '  gROOT->ProcessLine("pulls_compr()");' >> $compr_script
    echo '}' >> $compr_script
    root -l -n -q -b $compr_script >> output.txt
    rm $compr_script
    
done

