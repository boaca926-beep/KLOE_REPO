#!/bin/bash

compr_script=compr_script.C
echo '#include <iostream>' > $compr_script
echo "void compr_script() {" >> $compr_script
echo '  gROOT->ProcessLine(".L mass_compr_new.C");' >> $compr_script
echo "  gROOT->ProcessLine(mass_compr_new());" >> $compr_script

echo '  gROOT->ProcessLine(".L width_compr_new.C");' >> $compr_script
echo "  gROOT->ProcessLine(width_compr_new());" >> $compr_script

echo '  gROOT->ProcessLine(".L branch_compr_new.C");' >> $compr_script
echo "  gROOT->ProcessLine(branch_compr_new());" >> $compr_script
echo '}' >> $compr_script
root -l -n -q -b $compr_script
rm $compr_script
