# Run diff on key files in directory $1 to $2

if [[ -z $1 || -z $2 ]] ; then
    echo "Missing argument 1 and/or 2"
    echo "Usage: source diffKeyFiles.sh <srcDir> <diffDir>"
    echo "Example: source code/shell/diffKeyFiles.sh outTestHP11 goodTest"
    return
fi

mkdir -p $2

function checkIfFilesExistAndDiff {
    SRCDIR=$1
    DIFDIR=$2
    PATTERN=$3
    FILES=$SRCDIR/$PATTERN
    for f in $FILES
    do
        if [ ! -f $f ]; then
            echo "File $f does not exist, not diffed"
        else
            echo "Diffing $f"
            diff $f $DIFDIR/
        fi
    done
    }

checkIfFilesExistAndDiff $1 $2 'LCRs_*.tsv'
checkIfFilesExistAndDiff $1 $2 'BadKmers_*.tsv'
checkIfFilesExistAndDiff $1 $2 'IndelGroupsOverlapping_*.tsv'
checkIfFilesExistAndDiff $1 $2 'IndelGroupsNonoverlapping_*.tsv'
checkIfFilesExistAndDiff $1 $2 'NonvalidatedMarkers_*.tsv'
checkIfFilesExistAndDiff $1 $2 'MarkerErrors_*.tsv'
checkIfFilesExistAndDiff $1 $2 'MarkersOverlapping_*.tsv'
checkIfFilesExistAndDiff $1 $2 'MarkersNonoverlapping_*.tsv'
checkIfFilesExistAndDiff $1 $2 'MarkersAnnotated_*'
