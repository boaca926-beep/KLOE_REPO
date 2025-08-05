
#!/bin/bash

input=compr.h

FILE_TYPE=("vmd" "norm")
MODEL_TYPE=("VMD" "BW")

outtxt=fit_result.txt
echo "Fit results" > $outtxt

for (( j=0; j<1; j ++))
    
do
    echo "${FILE_TYPE[$j]}"

    file_type=\"${FILE_TYPE[$j]}\"
    model_type=\"${MODEL_TYPE[$j]}\"
    
    echo $file_type

    sed -i 's/\(const TString file_type =\)\(.*\)/\1 '$file_type';/' $input
    sed -i 's/\(const TString model_type =\)\(.*\)/\1 '$model_type';/' $input
    
    fit_script=fit_script.C
    echo '#include <iostream>' > $fit_script
    echo "void fit_script() {" >> $fit_script
    echo '  gROOT->ProcessLine(".L fit_result.C");' >> $fit_script
    #echo "  gROOT->ProcessLine(fit_result('${FILE_TYPE[$j]}'));" >> $fit_script
    echo "  gROOT->ProcessLine(fit_result());" >> $fit_script
    echo '}' >> $fit_script
    root -l -n -q -b $fit_script >> $outtxt

    plot_script=plot_script.C
    echo '#include <iostream>' > $plot_script
    echo "void plot_script() {" >> $plot_script
    echo '  gROOT->ProcessLine(".L plot_fit_results.C");' >> $plot_script
    echo "  gROOT->ProcessLine(plot_fit_results());" >> $plot_script
    echo '}' >> $plot_script
    root -l -n -q -b $plot_script  
    
done

rm $fit_script
rm $plot_script
    
