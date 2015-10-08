# $1: target for CLEANING
# Example: source code/shell/makeTestCase.sh ALL

if [[ -z $1 ]] ; then
    echo "Missing argument"
    echo "Usage: source makeTestCase.sh <CLEAN target>"
    echo "Example: source code/shell/makeTestCase.sh ALL"
    return
fi

make PARAMS=allParameters.test CLEAN=1 $1
make PARAMS=allParameters.test ALL | tee logFiles/makeLog_Test_$1.txt
make PARAMS=allParameters.test plotIndels | tee logFiles/makeLog_Test_Indels.txt
Rscript code/R/annotateFile.R annotate.template
Rscript code/R/annotateFile.R annotate/HP11_isInNearColumn.markers
Rscript code/R/annotateFile.R annotate/HP11_to_gff3.markers
