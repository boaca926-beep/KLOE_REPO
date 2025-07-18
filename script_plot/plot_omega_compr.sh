FILE_TYPE=("efficy_norm" "norm" "models" "Efficiency corrected" "Nominal")
#FILE_TYPE=("vertex" "norm" "omega" "" "")
#FILE_TYPE=("vmd_mini" "mini" "models")
#FILE_TYPE=("vmd_norm" "norm" "models" "VMD" "BW")
#FILE_TYPE=("isrlumi_norm" "norm" "models" "Exact" "Approx")

#sample_size="mini"
file_type1=${FILE_TYPE[0]}
file_type2=${FILE_TYPE[1]}
plot_type=${FILE_TYPE[2]}
legend_type1=${FILE_TYPE[3]}
legend_type2=${FILE_TYPE[4]}

# Plot main folder
plot_folder="../../plot_${plot_type}_compr/"
if [[ -d "$plot_folder" ]]; then
    ls $plot_folder
    echo "Remove plot folder"
    rm -rf $plot_folder
else
    echo "Create ${plot_folder}"
fi
mkdir $plot_folder

log_tmp="${plot_folder}log.txt"
#echo "" > ${log_tmp}
touch $log_tmp

header_compr=../header/plot_${plot_type}_compr.h

#echo $file_type1 $file_type2
echo $plot_type $legend_type1 $legend_type2 

sed -i 's|\(const TString file_type1 =\)\(.*\)|\1 "'"${file_type1}"'";|' "$header_compr"
sed -i 's|\(const TString file_type2 =\)\(.*\)|\1 "'"${file_type2}"'";|' "$header_compr"
sed -i 's|\(const TString outputFolder =\)\(.*\)|\1 "'"${plot_folder}"'";|' "$header_compr"

sed -i 's|\(const TString legend_type1 =\)\(.*\)|\1 "'"${legend_type1}"'";|' "$header_compr"
sed -i 's|\(const TString legend_type2 =\)\(.*\)|\1 "'"${legend_type2}"'";|' "$header_compr"

plot_compr_script=plot_compr_script.C
echo '#include <iostream>' > $plot_compr_script
echo "void plot_compr_script() {" >> $plot_compr_script
echo '  gROOT->ProcessLine(".L ../run_plot/plot_'${plot_type}'_compr.C");' >> $plot_compr_script
echo '  gROOT->ProcessLine("plot_'${plot_type}'_compr()");' >> $plot_compr_script
echo '}' >> $plot_compr_script
root -l -n -q -b $plot_compr_script >> ${log_tmp}

rm $plot_compr_script
