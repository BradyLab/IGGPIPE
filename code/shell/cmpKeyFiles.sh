# Run cmp on key files in directory $1 to $2

if [[ -z $1 || -z $2 ]] ; then
    echo "Missing argument 1 and/or 2"
    echo "Usage: source cmpKeyFiles.sh <srcDir> <cmpDir>"
    echo "Example: source code/shell/cmpKeyFiles.sh outTestHP11 goodTest"
    return
fi

mkdir -p $2

function checkIfFilesExistAndCmp {
    SRCDIR=$1
    CMPDIR=$2
    PATTERN=$3
    FILES=$SRCDIR/$PATTERN
    for f in $FILES
    do
        if [ ! -f $f ]; then
            echo "File $f does not exist, not compared"
        else
            echo "Comparing $f"
            BASENAME=$(basename "$f")
            cmp $f $CMPDIR/$BASENAME
        fi
    done
    }

checkIfFilesExistAndCmp $1 $2 'LCRs_*.tsv'
checkIfFilesExistAndCmp $1 $2 'BadKmers_*.tsv'
checkIfFilesExistAndCmp $1 $2 'IndelGroupsOverlapping_*.tsv'
checkIfFilesExistAndCmp $1 $2 'IndelGroupsNonoverlapping_*.tsv'
checkIfFilesExistAndCmp $1 $2 'NonvalidatedMarkers_*.tsv'
checkIfFilesExistAndCmp $1 $2 'MarkerErrors_*.tsv'
checkIfFilesExistAndCmp $1 $2 'MarkersOverlapping_*.tsv'
checkIfFilesExistAndCmp $1 $2 'MarkersNonoverlapping_*.tsv'
checkIfFilesExistAndCmp $1 $2 'MarkersAnnotated_*'
