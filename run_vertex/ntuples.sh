

ntuples_script=ntuples_script.C
    echo '#include <iostream>' > $ntuples_script
    echo "void ntuples_script() {" >> $ntuples_script
    echo 'gROOT->ProcessLine(".L ../run_vertex/ntuples.C");' >> $ntuples_script
    echo 'gROOT->ProcessLine("ntuples()");' >> $ntuples_script
    echo '}' >> $ntuples_script
    root -l -n -q -b $ntuples_script >> ${log_cut}
done
