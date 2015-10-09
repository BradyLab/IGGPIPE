# $1: target for make
# $2: target for make CLEAN=1  (default is same as $1, use 'NONE' to not clean)
# Example: source code/shell/makeTestCase.sh ALL

MAKETGT=$1
CLEANTGT=$2
if [[ -z $CLEANTGT ]] ; then
    CLEANTGT=$MAKETGT
fi

if [[ -z $MAKETGT ]] ; then
    echo "Usage: source makeTestCase.sh <make target> [<make CLEAN target>]"
    echo "Default for <make CLEAN target> is same as <make target>"
    echo "Use 'NONE' for <make CLEAN target> to not clean anything"
    echo "Example: source code/shell/makeTestCase.sh ALL"
    echo "Example: source code/shell/makeTestCase.sh ALL NONE"
    echo "Example: source code/shell/makeTestCase.sh findIndelGroups"
    return
fi

if [[ $CLEANTGT != "NONE" ]] ; then
    echo "make PARAMS=allParameters.test CLEAN=1 $CLEANTGT"
    make PARAMS=allParameters.test CLEAN=1 $CLEANTGT
fi

if [[ $CLEANTGT == "ALL" ]] ; then
    echo "make PARAMS=allParameters.test CLEAN=1 plotIndels"
    make PARAMS=allParameters.test CLEAN=1 plotIndels
    echo "rm -f $CLEANTGT/MarkersAnnotated_*"
    rm -f $CLEANTGT/MarkersAnnotated_*
fi

set -x #echo on

make PARAMS=allParameters.test $MAKETGT | tee logFiles/makeLog_Test_$MAKETGT.txt

if [[ $MAKETGT == "ALL" ]] ; then
    make PARAMS=allParameters.test plotIndels | tee logFiles/makeLog_Test_Indels.txt
    Rscript code/R/dotplot.R dotplot.template | tee logFiles/makeLog_Test_dotplot.txt
    Rscript code/R/annotateFile.R annotate.template | tee logFiles/makeLog_Test_annotate1.txt
    Rscript code/R/annotateFile.R annotate/HP11_isInNearColumn.markers | tee logFiles/makeLog_Test_annotate2.txt
    Rscript code/R/annotateFile.R annotate/HP11_to_gff3.markers | tee logFiles/makeLog_Test_annotate3.txt
fi

set +x #echo off
